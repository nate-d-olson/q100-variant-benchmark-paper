# Snakemake Workflow

This directory contains the Snakemake pipeline for processing GIAB v5q benchmark set files.

## Directory Structure

```
workflow/
├── Snakefile           # Main workflow entry point
├── rules/              # Rule definitions
│   ├── common.smk      # Helper functions and shared utilities
│   ├── downloads.smk   # Download and validation rules
│   ├── vcf_processing.smk  # VCF indexing and normalization
│   ├── var_tables.smk  # Variant annotation and table generation
│   └── exclusions.smk  # Exclusion region analysis
├── envs/               # Conda environment definitions
│   ├── bcftools.yaml
│   ├── bedtools.yaml
│   ├── downloads.yaml
│   ├── python.yaml
│   ├── rtg-tools.yaml
│   ├── samtools.yaml
│   └── truvari.yaml
└── scripts/            # Custom Python scripts
    ├── combine_beds_with_id.py
    ├── count_variants_by_type.py
    ├── expand_annotations.py
    ├── generate_header_lines.py
    └── generate_variant_parquet.py
```

## Rule Modules

All rule files are located in `workflow/rules/`. The pipeline includes 11 modular rule files:

1. **common.smk** (391 lines) - Helper functions and shared utilities
2. **downloads.smk** (329 lines) - Download and validation rules
3. **var_tables.smk** (233 lines) - Variant annotation and table generation
4. **exclusions.smk** (200 lines) - Exclusion region analysis
5. **strat_metrics.smk** (201 lines) - Stratification coverage metrics
6. **benchmark_comparisons.smk** (133 lines) - Benchmark set comparisons
7. **var_counts.smk** (110 lines) - Variant counting by type/stratification
8. **validation.smk** (87 lines) - Data validation (WIP - not yet integrated)
9. **vcf_processing.smk** (80 lines) - VCF indexing and normalization
10. **ref_processing.smk** (31 lines) - Reference genome preparation
11. **diff_tables.smk** (15 lines) - Regional coverage differences

### common.smk

Helper functions used across the pipeline:

- `get_var_table_inputs()` - Generate list of variant table output files
- `get_exclusion_*()` - Functions for exclusion file handling
- `get_stratification_*()` - Functions for stratification BED file handling
- `get_region_beds()` - Get benchmark region and exclusion BED paths
- `get_reference_checksum()` - Get checksum for reference genome validation
- `get_comparison_files()` - Get file paths for benchmark comparisons
- `get_bench_ids()`, `get_ref_ids()` - Get configured benchmark/reference IDs

### downloads.smk

Rules for downloading and validating remote input files:

| Rule                         | Description                                                 | Output                                                             |
| ---------------------------- | ----------------------------------------------------------- | ------------------------------------------------------------------ |
| `download_benchmark_vcf`     | Download benchmark VCF with SHA256 validation               | `resources/benchmarksets/{benchmark}_benchmark.vcf.gz`             |
| `download_benchmark_bed`     | Download benchmark BED with SHA256 validation               | `resources/benchmarksets/{benchmark}_benchmark.bed`                |
| `download_benchmark_dip_bed` | Download dipcall BED files                                  | `resources/benchmarksets/{benchmark}_dip.bed`                      |
| `prepare_reference`          | Download reference genome, convert to bgzip, create indices | `resources/references/{ref_name}.fa.gz`, `.fai`, `.gzi`            |
| `download_stratification`    | Download stratification BED files                           | `resources/stratifications/{ref}_{strat_name}.bed.gz`              |
| `download_exclusion`         | Download exclusion BED files                                | `resources/exclusions/{benchmark}/{exclusion_name}_{file_idx}.bed` |

### vcf_processing.smk

Rules for VCF processing and normalization:

| Rule                  | Description                                     | Output                                                 |
| --------------------- | ----------------------------------------------- | ------------------------------------------------------ |
| `index_vcf`           | Create tabix index for VCF files                | `{prefix}.vcf.gz.tbi`                                  |
| `split_multiallelics` | Split multiallelic variants using bcftools norm | `results/split_multiallelics/{benchmark}/split.vcf.gz` |

### var_tables.smk

Rules for variant annotation and table generation:

| Rule                           | Description                                | Output                                                                    |
| ------------------------------ | ------------------------------------------ | ------------------------------------------------------------------------- |
| `combine_stratification_beds`  | Merge stratification BEDs with IDs         | `results/combine_stratification_beds/{benchmark}/strat_combined.bed.gz`   |
| `combine_region_beds`          | Merge region BEDs (benchmark + exclusions) | `results/combine_region_beds/{benchmark}/region_combined.bed.gz`          |
| `run_truvari_anno_svinfo`      | Add SV annotations using truvari           | `results/run_truvari_anno_svinfo/{benchmark}/svinfo.vcf.gz`               |
| `generate_annotation_headers`  | Generate VCF header lines                  | `results/generate_annotation_headers/{benchmark}/annotation_headers.txt`  |
| `annotate_vcf_stratifications` | Annotate VCF with stratification IDs       | `results/annotate_vcf_stratifications/{benchmark}/strat_annotated.vcf.gz` |
| `annotate_vcf_regions`         | Annotate VCF with region IDs               | `results/annotate_vcf_regions/{benchmark}/fully_annotated.vcf.gz`         |
| `generate_variant_parquet`     | Generate Parquet table using Truvari       | `results/variant_tables/{benchmark}/variants.parquet`                     |

### exclusions.smk

Rules for exclusion region analysis:

| Rule                        | Description                                   | Output                                                             |
| --------------------------- | --------------------------------------------- | ------------------------------------------------------------------ |
| `materialize_exclusion`     | Process exclusion BED (merge pairs if needed) | `results/exclusions/{benchmark}/exclusions/{exclusion}.bed`        |
| `compute_dip_size`          | Calculate total benchmark region size         | `results/exclusions/{benchmark}/dip_size.txt`                      |
| `compute_exclusion_metrics` | Calculate overlap metrics for exclusions      | `results/exclusions/{benchmark}/metrics/{exclusion}.tsv`           |
| `aggregate_exclusion_table` | Combine all exclusion metrics into CSV        | `results/exclusions/{benchmark}/exclusions_intersection_table.csv` |

## Python Scripts

The pipeline includes 14 Python scripts in `workflow/scripts/`:

### Core Infrastructure

| Script              | Description                                                       |
| ------------------- | ----------------------------------------------------------------- |
| `logging_config.py` | Centralized logging configuration with structured output          |
| `exceptions.py`     | Custom exception classes (ValidationError, DataFormatError, etc.) |
| `validators.py`     | VCF/BED/TSV format validation utilities                           |

### Data Validation

| Script            | Description                              |
| ----------------- | ---------------------------------------- |
| `validate_vcf.py` | Validate VCF file format and structure   |
| `validate_bed.py` | Validate BED file format and coordinates |

### BED File Processing

| Script                    | Description                                        |
| ------------------------- | -------------------------------------------------- |
| `combine_beds_with_id.py` | Combine multiple BED files with unique identifiers |

### VCF Annotation & Querying

| Script                     | Description                                       |
| -------------------------- | ------------------------------------------------- |
| `generate_header_lines.py` | Generate VCF header lines for custom annotations  |
| `expand_annotations.py`    | Expand annotation ID lists to binary flag columns |

### Variant Analysis

| Script                       | Description                                         |
| ---------------------------- | --------------------------------------------------- |
| `generate_variant_parquet.py` | Generate Parquet variant tables using Truvari       |
| `count_variants_by_type.py`  | Count variants by type (SNV, INDEL, SV) from VCF   |
| `count_variants_by_strat.py` | Count variants per stratification region         |
| `stratify_comparison.py`     | Compare variants across stratifications          |

### Metrics Aggregation

| Script                      | Description                                  |
| --------------------------- | -------------------------------------------- |
| `summarize_var_counts.py`   | Aggregate variant count summaries            |
| `combine_metrics_counts.py` | Combine coverage metrics with variant counts |

## Conda Environments

| Environment      | Primary Tools                          |
| ---------------- | -------------------------------------- |
| `bcftools.yaml`  | bcftools for VCF manipulation          |
| `bedtools.yaml`  | bedtools for BED file operations       |
| `downloads.yaml` | wget for file downloads                |
| `python.yaml`    | Python with pysam for script execution |
| `rtg-tools.yaml` | RTG Tools for VCF statistics           |
| `samtools.yaml`  | samtools for reference genome indexing |
| `truvari.yaml`   | truvari for SV annotation              |

## Running the Pipeline

```bash
# Dry-run to validate workflow
snakemake -n --quiet

# Run with conda environments
snakemake --cores 4 --sdm conda

# Run specific target
snakemake results/variant_tables/v5q_grch38_smvar/variants.parquet --cores 4 --sdm conda
```

## Configuration

The pipeline is configured via `config/config.yaml` which defines:

- **benchmarksets**: Benchmark VCF/BED URLs and checksums
- **references**: Reference genome URLs, checksums, and stratification definitions
- **outputs**: Output directory paths

See `config/schema/config.schema.yaml` for the configuration schema.
