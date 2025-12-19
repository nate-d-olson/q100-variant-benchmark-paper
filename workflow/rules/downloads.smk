# Rules for downloading and verifying input files

import os

# Dynamically generate download rules for each reference defined in config
for ref_id, ref_info in config.get("references", {}).items():
    if "url" in ref_info and "path" in ref_info:

        rule:
            name:
                f"download_{ref_id}"
            output:
                ensure(ref_info["path"], **get_checksum_arg(ref_info.get("checksum"))),
            params:
                url=ref_info["url"],
            log:
                f"logs/downloads/{ref_id}.log",
            retries: 3
            conda:
                "../envs/bcftools.yaml"
            shell:
                "curl -L -o {output} {params.url} > {log} 2>&1"


rule prepare_reference:
    """Ensure reference is bgzipped and indexed for bcftools."""
    input:
        ref=lambda w: config["references"][w.ref]["path"],
    output:
        ref="resources/references/{ref}.fa.gz",
        fai="resources/references/{ref}.fa.gz.fai",
        gzi="resources/references/{ref}.fa.gz.gzi",
    log:
        "logs/references/{ref}_prepare.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        echo "Preparing reference {wildcards.ref}" > {log}

        # Check if input is bgzipped or gzipped
        if gzip -l {input.ref} > /dev/null 2>&1; then
            echo "Input is gzipped, decompressing and re-compressing with bgzip" >> {log}
            gunzip -c {input.ref} | bgzip -c > {output.ref} 2>> {log}
        else
            echo "Input is not compressed, compressing with bgzip" >> {log}
            bgzip -c {input.ref} > {output.ref} 2>> {log}
        fi

        # Index the reference
        samtools faidx {output.ref} 2>> {log}

        echo "Completed at $(date)" >> {log}
        """

# Dynamically generate download rules for benchmark sets
for benchmark, info in config.get("benchmarksets", {}).items():
    # Handle VCF download
    vcf_info = info.get("vcf")
    if isinstance(vcf_info, dict) and "url" in vcf_info:
        vcf_path = os.path.join(
            "resources/benchmarksets", get_filename_from_url(vcf_info["url"])
        )

        rule:
            name:
                f"download_benchmark_vcf_{benchmark}"
            output:
                ensure(vcf_path, **get_checksum_arg(vcf_info.get("checksum"))),
            params:
                url=vcf_info["url"],
            log:
                f"logs/downloads/benchmarks/{benchmark}_vcf.log",
            retries: 3
            conda:
                "../envs/bcftools.yaml"
            shell:
                "curl -L -o {output} {params.url} > {log} 2>&1"

    # Handle BED download
    bed_info = info.get("bed")
    if isinstance(bed_info, dict) and "url" in bed_info:
        bed_path = os.path.join(
            "resources/benchmarksets", get_filename_from_url(bed_info["url"])
        )

        rule:
            name:
                f"download_benchmark_bed_{benchmark}"
            output:
                ensure(bed_path, **get_checksum_arg(bed_info.get("checksum"))),
            params:
                url=bed_info["url"],
            log:
                f"logs/downloads/benchmarks/{benchmark}_bed.log",
            retries: 3
            conda:
                "../envs/bcftools.yaml"
            shell:
                "curl -L -o {output} {params.url} > {log} 2>&1"
