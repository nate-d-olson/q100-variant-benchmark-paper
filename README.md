# GIAB v5q HG002 Variant Benchmark Set Analysis

## Overview

This repository contains the analysis pipeline and manuscript for the **GIAB v5q HG002 Variant Benchmark Set**. The benchmark set provides high-confidence variant calls for HG002 (NA24385) using the T2T HG002 Q100 diploid assembly, with coordinates available for GRCh37, GRCh38, and CHM13v2.0 (T2T) reference genomes.

## Pipeline Purpose

The Snakemake pipeline processes benchmark set files to:

1. **Download and validate** benchmark VCF/BED files and reference genomes
2. **Annotate variants** with benchmark regions, exclusions (v5 only), and difficult genomic contexts using bcftools
3. **Generate variant tables** with comprehensive per-variant annotations (type, size, genomic context)
4. **Calculate overlap** with difficult genomic context (tandem repeats, homopolymers, segmental duplications, low mappability)
5. **Generate exclusion metrics** showing regions excluded from the benchmark and their impact
6. **Compare with historical benchmarks** (v4.2.1, v0.6 SV) to quantify improvements in callable regions and variant detection

## Key Outputs

- Annotated VCFs with region annotations (`results/annotate_vcf_regions/`)
- Comprehensive variant tables with per-variant annotations (`results/variant_tables/`)
- Exclusion intersection tables (`results/exclusions/`)
- **Benchmark Comparisons** (`results/stats/`):
  - `*_variants.csv`: Stratified counts of shared and unique variants
  - `*_regions.csv`: Stratified counts of shared and unique callable bases

## Analysis Resources

### Data Loading Functions

The `R/data_loading.R` module provides standardized functions for loading pipeline outputs:

- `load_genomic_context_metrics()` - Load primary analysis data with variant counts per genomic context (recommended)
- `load_primary_analysis_data()` - Load commonly used validated analysis data frames in one call
- `load_exclusion_metrics()` - Load exclusion overlaps (v5.0q only)
- `load_reference_sizes()` - Load reference genome sizes and N-content
- `load_variant_table()` - Load full variant data (use sparingly, files are large)
- `load_genomic_context_coverage()` - Load base-level coverage data
- `load_benchmark_regions()` - Load benchmark region BED files
- `load_platinum_pedigree_regions()` - Download and load Platinum Pedigree benchmark region BED files from public S3
- `parse_benchmark_id()` - Parse benchmark metadata from file path

All functions handle automatic file discovery, metadata parsing, and optional filtering.
Validation is performed during data loading before cache generation for cache-enabled
loaders, and schema validation is applied again when writing cache artifacts.

### Output Documentation

- **[Pipeline Outputs Reference](docs/pipeline-outputs.md)** - Complete catalog of output files with column schemas and examples
- **[Data Dictionary](docs/data-dictionary.md)** - Detailed definitions of metrics, stratifications, and variant types
- **[Output Relationships Diagram](docs/diagrams/output-relationships.mmd)** - Visual guide showing how files link together

# Getting Started

## Prerequisites

- [Mamba](https://mamba.readthedocs.io/) or [Conda](https://docs.conda.io/)
- [Snakemake](https://snakemake.readthedocs.io/) >= 8.0
- [snakefmt](https://github.com/snakemake/snakefmt) (for development)

## Installation

```bash
# Clone the repository
git clone https://github.com/nate-d-olson/q100-variant-benchmark-paper.git
cd q100-variant-benchmark-paper

# Create and activate the environment
mamba env create -f environment.yaml
mamba activate q100-smk
```

## Running the Pipeline

```bash
# Dry-run to validate workflow
snakemake -n --quiet

# Run pipeline locally with conda environments
snakemake --cores 4 --sdm conda

# Run using the Makefile shortcuts
make dry-run   # Validate workflow
make lint      # Check best practices
make test      # Run all validation tests
make run       # Execute the pipeline
```

## Testing with Local Fixtures

For fast, deterministic testing without downloading full datasets, the pipeline supports local file paths in configuration:

```bash
# Generate test fixtures (small subset of real data)
python scripts/create_grch38_debug_subset.py --force

# Run pipeline with test config using local fixtures
snakemake --configfile config/config.test_grch38_debug.yaml \
  --cores 4 --sdm conda \
  -- results/variant_tables/v5q_GRCh38_smvar/variants.parquet
```

**Local file support:**
- File paths with `file://` URLs (e.g., `file://$PWD/tests/fixtures/file.vcf.gz`)
- Absolute paths (e.g., `/path/to/file.vcf.gz`)
- Relative paths (e.g., `tests/fixtures/file.bed`)
- Environment variable expansion (e.g., `$PWD`, `$HOME`)
- All checksums still validated (SHA256/MD5)

See `config/config.test_grch38_debug.yaml` for a complete example using local fixtures.

## Repository Structure

```
project-root/
├── workflow/                        # Snakemake pipeline
│   ├── Snakefile                   # Main workflow entry point
│   ├── rules/                      # Rule definitions (13 files, 2,080 lines)
│   │   ├── common.smk              # Helper functions (523 lines)
│   │   ├── downloads.smk           # Download rules (329 lines)
│   │   ├── var_tables.smk          # Variant annotation & table generation (237 lines)
│   │   ├── genomic_context_metrics.smk  # Genomic context overlap metrics (155 lines)
│   │   ├── exclusions.smk          # Exclusion analysis (166 lines)
│   │   ├── benchmark_comparisons.smk    # Benchmark comparisons (137 lines)
│   │   ├── stratify_bench.smk      # Truvari stratification analysis (112 lines)
│   │   ├── var_counts.smk          # Variant counting (107 lines)
│   │   ├── output_organization.smk # Organize outputs (102 lines)
│   │   ├── validation.smk          # Data validation (87 lines)
│   │   ├── vcf_processing.smk      # VCF processing (79 lines)
│   │   ├── genomic_context_tables.smk   # Aggregate metrics (15 lines)
│   │   └── ref_processing.smk      # Reference indexing (31 lines)
│   ├── envs/                       # Conda environment definitions (8 files)
│   │   ├── bcftools.yaml
│   │   ├── bedtools.yaml
│   │   ├── downloads.yaml
│   │   ├── python.yaml
│   │   ├── rtg-tools.yaml
│   │   ├── samtools.yaml
│   │   ├── truvari.yaml
│   │   └── analysis.yaml
│   └── scripts/                    # Custom Python scripts (15 files)
│       ├── logging_config.py
│       ├── exceptions.py
│       ├── validators.py
│       ├── validate_vcf.py
│       ├── validate_bed.py
│       ├── combine_beds_with_id.py
│       ├── compute_bed_metrics.py
│       ├── extract_info_fields.py
│       ├── generate_header_lines.py
│       ├── count_variants_by_genomic_context.py
│       ├── aggregate_stratified_bench.py
│       ├── stratify_comparison.py
│       ├── summarize_var_counts.py
│       └── combine_metrics_counts.py
├── config/
│   ├── config.yaml                 # Pipeline configuration
│   └── schema/
│       └── config.schema.yaml      # Configuration validation schema
├── docs/                           # Comprehensive documentation
│   ├── architecture.md             # System architecture
│   ├── api-reference.md            # Python API reference
│   └── troubleshooting.md          # Common issues and solutions
├── tests/                          # Unit tests
│   ├── unit/
│   │   ├── test_validators.py
│   │   └── test_exceptions.py
│   └── fixtures/
│       ├── sample.vcf
│       └── sample.bed
├── manuscript/                     # Quarto manuscript files
│   ├── introduction.qmd
│   ├── methods.qmd
│   ├── results.qmd
│   ├── discussion.qmd
│   └── references.bib
├── analysis/                       # Supplementary analysis notebooks
├── resources/                      # Downloaded resources (gitignored)
│   ├── benchmarksets/
│   ├── references/
│   ├── stratifications/
│   └── exclusions/
├── results/                        # Pipeline outputs (gitignored)
│   ├── variant_tables/
│   ├── var_counts/
│   ├── genomic_context_metrics/
│   ├── exclusions/
│   └── stats/
├── environment.yaml                # Mamba environment spec
├── pyproject.toml                  # Python project configuration
├── Makefile                        # Development shortcuts
└── _quarto.yml                     # Quarto project config
```

# Benchmark Sets Analyzed

| Benchmark ID          | Reference | Variant Type        | Description                              |
| --------------------- | --------- | ------------------- | ---------------------------------------- |
| v5.0q_CHM13v2.0_smvar | CHM13v2.0 | Small variants      | T2T Q100 benchmark (SNVs + indels <50bp) |
| v5.0q_CHM13v2.0_stvar | CHM13v2.0 | Structural variants | T2T Q100 benchmark (SVs ≥50bp)           |
| v5.0q_GRCh38_smvar    | GRCh38    | Small variants      | T2T Q100 liftover                        |
| v5.0q_GRCh38_stvar    | GRCh38    | Structural variants | T2T Q100 liftover                        |
| v5.0q_GRCh37_smvar    | GRCh37    | Small variants      | T2T Q100 liftover                        |
| v5.0q_GRCh37_stvar    | GRCh37    | Structural variants | T2T Q100 liftover                        |
| v4.2.1_GRCh38_smvar   | GRCh38    | Small variants      | Historical GIAB v4.2.1                   |
| v0.6_GRCh37_stvar     | GRCh37    | Structural variants | Historical SV v0.6                       |

# Genomic Stratification Contexts

The pipeline analyzes benchmark overlap with GIAB v3.6 stratification regions:

| Context | Description                                      |
| ------- | ------------------------------------------------ |
| TR      | All tandem repeats                               |
| TR10kb  | Tandem repeats ≥10001bp with 5bp slop            |
| HP      | Homopolymers ≥7bp (perfect) or ≥11bp (imperfect) |
| SD      | All segmental duplications                       |
| SD10kb  | Segmental duplications >10kb                     |
| MAP     | Low mappability regions                          |

# Pipeline Improvements (2026-01)

Recent enhancements have significantly improved pipeline reliability, maintainability, and usability. See [`IMPROVEMENT_SUGGESTIONS.md`](IMPROVEMENT_SUGGESTIONS.md) and [`IMPLEMENTATION_SUMMARY.md`](IMPLEMENTATION_SUMMARY.md) for complete details.

### Comprehensive Documentation

- **[Architecture Documentation](docs/architecture.md)** - System design with visual diagrams
- **[API Reference](docs/api-reference.md)** - Complete documentation of 15+ helper functions
- **[Troubleshooting Guide](docs/troubleshooting.md)** - Solutions for common issues

## Development

### Code Quality

```bash
# Run linting and format checks
make test

# Format Snakemake files
make format
```

### GitHub Copilot Support

This repository is configured for optimal use with GitHub Copilot Coding Agent:

- **Comprehensive Instructions**: Custom instructions in `.github/copilot-instructions.md` and `.github/AGENTS.md`
- **Path-Specific Guidance**: Specialized instructions for Snakemake rules, Python scripts, and R analysis code
- **Automated Setup**: Environment setup workflow in `.github/workflows/copilot-setup-steps.yml`

For details, see:
- **[Copilot Setup Verification](.github/COPILOT_SETUP_VERIFICATION.md)** - Complete setup documentation
- **[Instructions Directory](.github/instructions/)** - Path-specific coding guidelines

### Commit Conventions

This project uses [Conventional Commits](https://www.conventionalcommits.org/):

- `feat:` - New feature or analysis
- `fix:` - Bug fix
- `docs:` - Documentation updates
- `chore:` - Maintenance (dependencies, configs)
- `refactor:` - Code restructuring

## Citation

If you use this benchmark set, please cite:

> [Citation to be added upon publication]

# License

This project is licensed under CC0 1.0 Universal - see the LICENSE file for details.
