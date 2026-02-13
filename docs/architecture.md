# Q100 Variant Benchmark Pipeline Architecture

## Overview

The Q100 variant benchmark pipeline is a Snakemake-based workflow for analyzing and characterizing the GIAB v5q HG002 variant benchmark set. The pipeline processes high-confidence variant calls across multiple reference genomes and compares them with historical benchmarks.

## System Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                    Configuration Layer                          │
│  config/config.yaml - Benchmark definitions, references, URLs   │
└────────────────┬────────────────────────────────────────────────┘
                 │
                 ▼
┌─────────────────────────────────────────────────────────────────┐
│                      Workflow Engine                            │
│              Snakemake 8.0+ (DAG execution)                     │
└────┬────────────────────────────────────────────────────────┬───┘
     │                                                        │
     ▼                                                        ▼
┌─────────────────┐                                  ┌────────────────┐
│  Download Phase │                                  │ Validation     │
│  - Benchmarks   │                                  │ Phase          │
│  - References   │                                  │ - VCF checks   │
│  - Stratifications│──────────────────────────────▶ │ - BED checks   │
│  - Exclusions   │                                  │ - Format verify│
└────┬────────────┘                                  └────────┬───────┘
     │                                                        │
     ▼                                                        │
┌─────────────────────────────────────────────────────────────┴───┐
│                   Processing Phase                              │
│  ┌───────────────┐  ┌────────────────┐  ┌────────────────────┐ │
│  │ VCF Processing│  │ BED Processing │  │ Reference Indexing │ │
│  │ - Index       │  │ - Combine      │  │ - FAI creation     │ │
│  │ - Normalize   │  │ - Annotate IDs │  │                    │ │
│  └───────────────┘  └────────────────┘  └────────────────────┘ │
└────────────────────────────┬────────────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────────────┐
│                    Annotation Phase                             │
│  ┌──────────────────────────────────────────────────────────┐   │
│  │ Combine stratification BEDs with unique IDs              │   │
│  │ Input: Multiple BED files per reference                 │   │
│  │ Output: Single BED with STRAT_ID column                 │   │
│  └──────────────────────────────────────────────────────────┘   │
│  ┌──────────────────────────────────────────────────────────┐   │
│  │ Combine exclusion/region BEDs with unique IDs           │   │
│  │ Input: Benchmark regions + exclusion categories         │   │
│  │ Output: Single BED with REGION_ID column                │   │
│  └──────────────────────────────────────────────────────────┘   │
└────────────────────────────┬────────────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────────────┐
│                      Query Phase                                │
│  ┌──────────────────────────────────────────────────────────┐   │
│  │ bcftools query - Extract variant INFO fields            │   │
│  │ - Chromosome, position, type                            │   │
│  │ - Annotations (STRAT_IDS, REGION_IDS)                   │   │
│  │ Output: TSV with one row per variant                    │   │
│  └──────────────────────────────────────────────────────────┘   │
└────────────────────────────┬────────────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────────────┐
│                  Table Generation Phase                         │
│  ┌──────────────────────────────────────────────────────────┐   │
│  │ Expand annotations - Convert ID lists to binary flags   │   │
│  │ Input: TSV with STRAT_IDS, REGION_IDS columns           │   │
│  │ Output: TSV with STRAT_* and REGION_* binary columns    │   │
│  └──────────────────────────────────────────────────────────┘   │
└────────────────────────────┬────────────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────────────┐
│                     Metrics Phase                               │
│  ┌────────────────────┐  ┌─────────────────────────────────┐   │
│  │ Coverage Metrics   │  │ Variant Counts                  │   │
│  │ - Per stratification│  │ - By type (SNV, INDEL, SV)     │   │
│  │ - Overlap stats    │  │ - By stratification             │   │
│  └────────────────────┘  └─────────────────────────────────┘   │
└────────────────────────────┬────────────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────────────┐
│                  Aggregation Phase                              │
│  - Combine metrics into summary tables                          │
│  - Generate comparison reports                                  │
│  - Create visualizations (via Quarto notebooks)                 │
└─────────────────────────────────────────────────────────────────┘
```

## Data Flow

### Input Sources

1. **Benchmark VCF Files**
   - GIAB v5q variants (SNVs, INDELs, SVs)
   - Historical benchmarks (v4.2.1, v0.6 SV)
   - Coordinates for GRCh37, GRCh38, CHM13v2.0

2. **Benchmark BED Files**
   - High-confidence regions for each benchmark

3. **Reference Genomes**
   - GRCh37, GRCh38, CHM13v2.0 FASTA files

4. **Stratification BEDs**
   - Tandem repeats
   - Homopolymers (7bp+)
   - Segmental duplications
   - Low mappability regions
   - Difficult-to-sequence regions
   - Other genomic features

5. **Exclusion BEDs**
   - Regions excluded from benchmarks
   - Multiple exclusion categories per benchmark

### Processing Pipeline

```
Benchmarks (VCF) ──┐
References (FASTA) ─┼─▶ [Validation] ─▶ [Index/Normalize] ─┐
Stratifications ────┤                                       │
Exclusions ─────────┘                                       ▼
                                                    [Annotate with IDs]
                                                            │
                                                            ▼
                                              [Query variant + annotations]
                                                            │
                                                            ▼
                                                [Expand to binary flags]
                                                            │
                                                            ▼
                                                    [Calculate metrics]
                                                            │
                                                            ▼
                                                    [Aggregate results]
```

### Output Structure

```
results/
├── validation/               # Validation reports
│   ├── {benchmark}/{ref}/
│   └── stratifications/{ref}/
├── benchmarks/               # Indexed VCFs and regions
│   └── {benchmark}/{ref}/
├── references/               # Indexed reference genomes
│   └── {ref}/
├── var_tables/               # Variant annotation tables
│   └── {benchmark}/{ref}/
├── metrics/                  # Stratification metrics
│   └── {benchmark}/{ref}/
├── var_counts/               # Variant count summaries (PRIMARY OUTPUTS)
│   └── {benchmark}/{ref}/
│       └── genomic_context_combined_metrics.csv ✓ Use first
├── exclusions/               # Exclusion overlaps (v5.0q only)
│   └── {benchmark}/{ref}/
│       └── exclusions_intersection_table.csv ✓ Use first
├── variant_tables/           # Full variant data (DETAILED)
│   └── {benchmark}/{ref}/
│       └── variants.tsv ⚠️ Large, use filters
├── diff_region_coverage/     # Base-level coverage (DETAILED)
│   └── {benchmark}/{ref}/
│       ├── HP_cov.bed ⚠️ Large
│       ├── MAP_cov.bed
│       ├── SD_cov.bed
│       ├── SD10kb_cov.bed
│       ├── TR_cov.bed
│       └── TR10kb_cov.bed
├── ref_genome_sizes/         # Reference metadata (PRIMARY)
│   ├── GRCh37_size.tsv ✓ Use first
│   ├── GRCh38_size.tsv
│   └── CHM13v2.0_size.tsv
└── aggregated/               # Final summary tables
```

**Legend:**

- ✓ = Primary outputs (small, fast to load) - recommended for most analyses
- ⚠️ = Detailed outputs (large, slow) - use only when variant-level detail needed

### Primary vs. Detailed Output Files

The pipeline generates two tiers of outputs optimized for different use cases:

#### 1. Primary Analysis Files (Small, ~KB per file)

- `genomic_context_combined_metrics.csv` - Aggregated metrics with variant counts
- `exclusions_intersection_table.csv` - Exclusion overlap summaries (v5.0q only)
- `ref_genome_sizes/*.tsv` - Reference genome metadata

**Characteristics:**

- Fast to load (seconds, not minutes)
- Pre-aggregated for common analyses
- Sufficient for 95% of analyses
- Generated by combining intermediate results

**When to use:**

- Creating comparison plots and summary tables
- Analyzing benchmark completeness
- Calculating variant densities and coverage
- All quick exploratory analyses

#### 2. Detailed Data Files (Large, ~MB-GB per file)

- `variant_tables/*/variants.tsv` - Full variant-level annotations
- `diff_region_coverage/*/*.bed` - Base-level coverage data

**Characteristics:**

- Slow to load (minutes, high memory)
- Variant-level or base-level granularity
- Required only for in-depth investigations
- Use with filters to reduce memory usage

**When to use:**

- Analyzing specific variant properties
- Investigating edge cases or anomalies
- Variant-by-variant comparisons
- Base-level coverage analysis

### Recommended Analysis Workflow

1. **Start with primary files** using `R/data_loading.R` functions
2. **Generate summary statistics** and visualizations
3. **If variant-level detail needed**, load detailed files with filters:

   ```r
   # Load specific chromosomes or variant types to reduce memory
   vars <- load_variant_table(
       "v5.0q_GRCh38_smvar",
       filters = list(
           chromosomes = c("chr1", "chr2"),
           variant_types = c("SNV"),
           in_benchmark_only = TRUE
       )
   )
   ```

4. **Document data sources** used for each analysis result

## Module Breakdown

### Core Workflow Rules

The pipeline consists of 13 modular rule files totaling 2,080 lines of Snakemake code:

| Module                        | Lines | Purpose                                  | Key Rules                                                                  |
| ----------------------------- | ----- | ---------------------------------------- | -------------------------------------------------------------------------- |
| `common.smk`                  | 523   | Helper functions and config parsing      | `get_region_beds()`, `get_exclusion_inputs()`, `get_genomic_context_ids()` |
| `downloads.smk`               | 329   | File downloads with SHA256 validation    | `download_benchmark_vcf`, `download_reference`                             |
| `var_tables.smk`              | 237   | Variant annotation and table generation  | `combine_genomic_context_beds`, `annotate_vcf_genomic_contexts`            |
| `genomic_context_metrics.smk` | 155   | Genomic context coverage metrics         | `compute_genomic_context_metrics`, `aggregate_genomic_context_metrics`     |
| `exclusions.smk`              | 166   | Exclusion region analysis                | `materialize_exclusion`, `compute_exclusion_metrics`                       |
| `benchmark_comparisons.smk`   | 137   | Benchmark set comparisons                | `run_truvari_compare`, `stratify_comparison`                               |
| `stratify_bench.smk`          | 112   | Truvari stratification analysis          | `truvari_stratify`                                                         |
| `var_counts.smk`              | 107   | Variant counting by type/genomic context | `count_variants_by_genomic_context`, `combine_metrics_and_counts`          |
| `output_organization.smk`     | 102   | Organize outputs into final structure    | `organize_outputs`                                                         |
| `validation.smk`              | 87    | Data validation                          | `validate_benchmark_vcf`, `validate_benchmark_bed`                         |
| `vcf_processing.smk`          | 79    | VCF indexing and normalization           | `index_vcf`, `split_multiallelics`                                         |
| `genomic_context_tables.smk`  | 15    | Aggregate genomic context tables         | `aggregate_genomic_context_tables`                                         |
| `ref_processing.smk`          | 31    | Reference genome indexing                | `index_reference`                                                          |

### Python Data Processing Scripts

The pipeline includes 15 Python scripts organized by function:

| Script                                 | Purpose                                | Key Features                           |
| -------------------------------------- | -------------------------------------- | -------------------------------------- |
| **Core Infrastructure**                |                                        |                                        |
| `logging_config.py`                    | Centralized logging setup              | Structured logs, multi-handler support |
| `exceptions.py`                        | Custom exception classes               | Context-rich error messages            |
| `validators.py`                        | Data validation utilities              | VCF/BED/TSV format checking            |
| **Data Validation**                    |                                        |                                        |
| `validate_vcf.py`                      | VCF format validation                  | Pre-flight format checks               |
| `validate_bed.py`                      | BED format validation                  | Coordinate and structure checks        |
| **BED Processing**                     |                                        |                                        |
| `combine_beds_with_id.py`              | Merge BED files with unique IDs        | Adds ID column for tracking            |
| `compute_bed_metrics.py`               | Compute coverage metrics from BEDs     | Intersect and percentage calculations  |
| **VCF Annotation**                     |                                        |                                        |
| `extract_info_fields.py`               | Extract VCF INFO field names           | Dynamic header generation              |
| `generate_header_lines.py`             | Create VCF annotation headers          | For bcftools annotate                  |
| **Variant Analysis**                   |                                        |                                        |
| `count_variants_by_genomic_context.py` | Count variants per genomic context     | Aggregates variant types               |
| `aggregate_stratified_bench.py`        | Aggregate stratified benchmark results | Cross-context summaries                |
| `stratify_comparison.py`               | Compare variants across benchmarks     | Benchmark comparison analysis          |
| **Metrics Aggregation**                |                                        |                                        |
| `summarize_var_counts.py`              | Aggregate variant count tables         | Cross-benchmark summaries              |
| `combine_metrics_counts.py`            | Combine metrics into final tables      | Merges coverage + counts               |

## Key Design Patterns

### 1. Configuration-Driven Design

Benchmarks are defined in `config/config.yaml` with:

- Unique benchmark IDs
- Reference genome mappings
- URLs and SHA256 checksums
- Stratification types per reference
- Exclusion categories

**Benefit**: Easy to add new benchmarks without code changes.

### 2. Modular Rule Organization

Rules are organized by function into separate `.smk` files:

- Clear separation of concerns
- Independent development and testing
- Easy to understand dependencies

### 3. SHA256 Validation

All downloads include checksum verification:

- Ensures data integrity
- Prevents corrupt file propagation
- Reproducible results

### 4. Conda Environment Isolation

Each rule type uses dedicated conda environments:

- `base.yaml` - Python utilities
- `bedtools.yaml` - BED operations
- `bcftools.yaml` - VCF operations
- `samtools.yaml` - Reference indexing

## Complex Logic Explanation

### `get_region_beds()` Function

This function (in `common.smk`) combines benchmark regions with exclusion categories:

```python
def get_region_beds(benchmark, ref):
    """
    Get all BED files for a benchmark's regions and exclusions.

    Returns:
        List of BED file paths with format:
        [
            (benchmark_regions.bed, "BMKREGIONS"),
            (exclusion1.bed, "EXCL_CATEGORY1"),
            (exclusion2.bed, "EXCL_CATEGORY2"),
            ...
        ]
    """
```

**Why it's complex**: Handles multiple exclusion categories with dynamic naming based on config.

**Solution**: String manipulation to convert exclusion names ("consecutive-svs") to IDs ("CONSECUTIVE_SVS").

### Annotation Expansion

The `expand_annotations.py` script converts comma-separated ID lists into binary columns:

**Input**:

```
CHROM  POS  STRAT_IDS          REGION_IDS
chr1   100  TANDEM_REPEATS,SD  BMKREGIONS,EXCL_LCR
```

**Output**:

```
CHROM  POS  STRAT_TANDEM_REPEATS  STRAT_SD  BMKREGIONS  EXCL_LCR
chr1   100  1                     1         1           1
```

**Why needed**: Binary flags enable efficient filtering and aggregation in downstream analysis.

## Dependency Graph

```
config.yaml
    │
    ├─▶ Download resources
    │       │
    │       ├─▶ Validate formats (NEW)
    │       │       │
    │       │       └─▶ Index VCFs/References
    │       │               │
    │       │               └─▶ Combine BEDs with IDs
    │       │                       │
    │       │                       └─▶ Annotate VCF
    │       │                               │
    │       │                               └─▶ Query variants
    │       │                                       │
    │       │                                       └─▶ Expand annotations
    │       │                                               │
    │       │                                               └─▶ Calculate metrics
    │       │                                                       │
    │       │                                                       └─▶ Aggregate results
    │       │
    └───────┴─▶ Analysis notebooks
                    │
                    └─▶ Manuscript figures
```

## R Data Loading and Caching Layer

The R analysis layer provides data loading, schema validation, and Parquet-based caching for pipeline outputs.

### Module Architecture

```
R/
├── schemas.R       ← Arrow schema registry, factor levels, validation rules
├── cache.R         ← Parquet caching infrastructure (sources schemas.R)
└── data_loading.R  ← Loading functions with cache integration (sources both)
```

### Schema Registry (`R/schemas.R`)

Centralized definitions used by both writing and reading:

| Component                       | Purpose                                             |
| ------------------------------- | --------------------------------------------------- |
| `get_arrow_schema(dataset)`     | Arrow type definitions for Parquet write            |
| `get_factor_levels(dataset)`    | Factor level restoration after Parquet read         |
| `get_validation_rules(dataset)` | Column predicate functions for fail-fast validation |
| `validate_data(df, dataset)`    | Run all validations on a data frame                 |

Supported datasets: `variant_table`, `diff_coverage`, `benchmark_regions`

### Caching Infrastructure (`R/cache.R`)

Parquet-only caching with automatic invalidation:

| Function                      | Purpose                                                |
| ----------------------------- | ------------------------------------------------------ |
| `write_cache()`               | Validate, strip factors, embed metadata, write Parquet |
| `read_cache()`                | Read Parquet, restore factors, return tibble (or NULL) |
| `cache_is_valid()`            | Check if valid cache exists for given inputs           |
| `invalidate_cache(dataset)`   | Remove all cache files for a dataset                   |
| `clear_cache()`               | Remove entire cache directory                          |
| `cache_info()`                | List cached files with size, age, dataset name         |
| `collect_pipeline_metadata()` | Gather R version, package versions, config summary     |

**Design decisions:**

- Cache directory configurable via `getOption("q100.cache_dir")` (default: `analysis/cache/`)
- Cache key: `rlang::hash()` of dataset name + source file mtimes + parameters
- Compression: zstd level 3 for speed/compression tradeoff
- Pipeline metadata stored as JSON in Parquet file-level key-value metadata
- Factors stripped before write, restored on read from schema registry
- Validation runs on write (fail-fast), optional on read

### Data Loading Functions (`R/data_loading.R`)

Three functions support caching via `use_cache` and `force_refresh` parameters:

| Function                          | Dataset             | Source Files            |
| --------------------------------- | ------------------- | ----------------------- |
| `load_variant_table()`            | `variant_table`     | `variants.tsv`          |
| `load_genomic_context_coverage()` | `diff_coverage`     | `*_cov.bed` files       |
| `load_benchmark_regions()`        | `benchmark_regions` | `*_benchmark.bed` files |

Functions NOT cached (small datasets, load directly from pipeline output):

- `load_genomic_context_metrics()`, `load_exclusion_metrics()`, `load_reference_sizes()`

### Adding a New Cached Dataset

1. **`R/schemas.R`**: Add entries to `get_arrow_schema()`, `get_factor_levels()`, `get_validation_rules()`
2. **`R/data_loading.R`**: Add `use_cache`/`force_refresh` params, call `read_cache()` before loading, `write_cache()` after

No changes needed in `R/cache.R` -- it is dataset-agnostic.

## Technology Stack

| Layer           | Technology          | Purpose                                             |
| --------------- | ------------------- | --------------------------------------------------- |
| Workflow        | Snakemake 8.0+      | DAG execution, rule management                      |
| Environment     | Conda/Mamba         | Reproducible package management                     |
| VCF             | bcftools            | VCF query, annotation, manipulation                 |
| BED             | bedtools            | Genomic region operations                           |
| Indexing        | samtools            | FASTA indexing                                      |
| Data Processing | Python 3.x          | Custom scripts with logging (NEW)                   |
| Analysis        | R 4.5 + Quarto      | Statistical analysis and reporting                  |
| Data Caching    | Arrow/Parquet       | Schema-validated Parquet caching for large datasets |
| Code Quality    | air + lintr         | R code formatting and linting                       |
| Validation      | Custom Python (NEW) | Format checking, integrity validation               |

## Error Handling (NEW)

The pipeline now includes comprehensive error handling:

1. **Structured Logging**
   - Consistent format: `[TIMESTAMP] [LEVEL] [MODULE] Message`
   - Captured by Snakemake in `logs/` directory

2. **Custom Exceptions**
   - `ValidationError` - File existence, format, integrity
   - `DataFormatError` - VCF/BED/TSV structure errors
   - `ProcessingError` - Data transformation failures
   - `ConfigurationError` - Invalid pipeline configuration

3. **Validation Layer**
   - Pre-flight checks before expensive operations
   - Early detection of corruption or format violations
   - Detailed validation reports in `results/validation/`

## Performance Characteristics

| Phase               | Bottleneck        | Typical Runtime | Memory       |
| ------------------- | ----------------- | --------------- | ------------ |
| Download            | Network bandwidth | Varies          | Low          |
| Validation (NEW)    | I/O               | 1-5 min         | Low          |
| VCF Indexing        | I/O               | 5-10 min        | Low          |
| BED Combination     | I/O               | 1-2 min         | Moderate     |
| VCF Annotation      | bcftools          | 10-30 min       | Moderate     |
| Metrics Calculation | bedtools          | 5-15 min        | Low-Moderate |

**Parallelization**: Pipeline supports `-j N` for parallel rule execution.

## Recent Improvements (This Branch)

1. **Structured Logging Framework**
   - `logging_config.py` - Centralized configuration
   - `exceptions.py` - Context-rich error messages
   - Updated 2 scripts with logging (8 total to update)

2. **Data Validation Layer**
   - `validators.py` - VCF/BED/TSV validation utilities
   - `validate_vcf.py`, `validate_bed.py` - Validation scripts
   - `validation.smk` - Snakemake validation rules

3. **Type Hints and Documentation**
   - Added type hints to updated Python scripts
   - Improved docstrings with examples
   - Better error messages with actionable suggestions

4. **Parquet Caching and Schema Registry**
   - `R/schemas.R` - Arrow schema definitions, factor levels, validation rules
   - `R/cache.R` - Parquet caching with zstd compression and pipeline metadata
   - Cache integration in `load_variant_table()`, `load_genomic_context_coverage()`, `load_benchmark_regions()`
   - 45 tests in `tests/test_cache.R` covering schemas, cache round-trips, invalidation, metadata

5. **Code Quality Tooling**
   - `air.toml` - Air formatter configuration (100-char lines, 2-space indent)
   - `.lintr` - Lintr configuration (snake_case + dotted.case, 100-char lines)

## Future Enhancements

See `IMPROVEMENT_SUGGESTIONS.md` for:

- Unit testing framework (pytest)
- Architecture documentation expansion (this file!)
- Performance profiling and optimization
- Additional Python script updates with logging

---

_Last Updated: 2026-02-07_
_Documentation synchronized with codebase (common.smk duplicate functions removed)_
