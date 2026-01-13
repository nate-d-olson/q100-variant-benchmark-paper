"""
Data validation rules for Q100 variant benchmark pipeline.

Validates input files and intermediate outputs to catch errors early.
"""


rule validate_benchmark_vcf:
    """Validate benchmark VCF file format and structure."""
    input:
        vcf="resources/benchmarks/{benchmark}/{ref}/variants.vcf.gz",
    output:
        report="results/validation/{benchmark}/{ref}/vcf_validation.txt",
    log:
        "logs/validation/{benchmark}/{ref}/vcf.log",
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/validate_vcf.py"


rule validate_benchmark_bed:
    """Validate benchmark BED file format and structure."""
    input:
        bed="resources/benchmarks/{benchmark}/{ref}/regions.bed",
    output:
        report="results/validation/{benchmark}/{ref}/bed_validation.txt",
    log:
        "logs/validation/{benchmark}/{ref}/bed.log",
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/validate_bed.py"


rule validate_stratification_beds:
    """Validate stratification BED files."""
    input:
        bed="resources/stratifications/{ref}/{strat_type}.bed.gz",
    output:
        report="results/validation/stratifications/{ref}/{strat_type}_validation.txt",
    log:
        "logs/validation/stratifications/{ref}/{strat_type}.log",
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/validate_bed.py"


rule validation_summary:
    """Generate summary report of all validation results."""
    input:
        # Aggregate all validation outputs
        vcf_reports=expand(
            "results/validation/{benchmark}/{ref}/vcf_validation.txt",
            benchmark=get_benchmarks(),
            ref=lambda wc: get_refs_for_benchmark(wc.benchmark)
        ),
        bed_reports=expand(
            "results/validation/{benchmark}/{ref}/bed_validation.txt",
            benchmark=get_benchmarks(),
            ref=lambda wc: get_refs_for_benchmark(wc.benchmark)
        ),
    output:
        summary="results/validation/validation_summary.txt",
    log:
        "logs/validation/summary.log",
    shell:
        """
        echo "Validation Summary - $(date)" > {output.summary}
        echo "=" >> {output.summary}
        echo "" >> {output.summary}

        echo "VCF Validations:" >> {output.summary}
        for report in {input.vcf_reports}; do
            echo "  - $report: $(tail -n1 $report)" >> {output.summary}
        done

        echo "" >> {output.summary}
        echo "BED Validations:" >> {output.summary}
        for report in {input.bed_reports}; do
            echo "  - $report: $(tail -n1 $report)" >> {output.summary}
        done

        echo "" >> {output.summary}
        echo "Validation complete!" >> {output.summary}
        """
