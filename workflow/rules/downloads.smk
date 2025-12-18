# Rules for downloading and verifying input files

import os


# Function to determine checksum type and return appropriate ensure() argument
def get_checksum_arg(checksum):
    if not checksum:
        return {}
    # MD5 is 32 hex chars, SHA256 is 64 hex chars
    if len(checksum) == 32:
        return {"md5": checksum}
    elif len(checksum) == 64:
        return {"sha256": checksum}
    else:
        # Default to md5 if length doesn't match standard sizes, or let Snakemake handle it
        return {"md5": checksum}


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
