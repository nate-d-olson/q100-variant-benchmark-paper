# Data Dictionary for R Loading Functions

This document describes the data objects returned by the functions in `R/data_loading.R`. It serves as a reference for development and analysis.

> **Note:** This document must be updated whenever changes are made to `R/data_loading.R` or `R/schemas.R` that affect the structure of the returned data objects.

## 1. Genomic Context Metrics

**Function:** `load_genomic_context_metrics()`

**Description:** Primary analysis data containing aggregated metrics and variant counts per genomic context (e.g., Homopolymers, Tandem Repeats). This is the most commonly used dataset for high-level comparisons.

| Column | Type | Description |
| :--- | :--- | :--- |
| `bench_version` | Factor | Benchmark version (`v0.6`, `v4.2.1`, `v5.0q`) |
| `ref` | Factor | Reference genome (`GRCh37`, `GRCh38`, `CHM13v2.0`) |
| `bench_type` | Factor | Benchmark set type (`smvar`, `stvar`) |
| `context_name` | Factor | Genomic context identifier (e.g., `HP`, `TR`, `SD`) |
| `context_bp` | Numeric | Total size of the genomic context in bases |
| `intersect_bp` | Numeric | Number of bases where the genomic context overlaps the benchmark regions |
| `pct_of_context` | Numeric | Percentage of the genomic context covered by the benchmark |
| `pct_of_bench` | Numeric | Percentage of the benchmark covered by this genomic context |
| `total_variants` | Integer | Total number of variants in this context |
| `snv_count` | Integer | Count of Single Nucleotide Variants |
| `indel_count` | Integer | Count of Insertions/Deletions |
| `del_count` | Integer | Count of Deletions (Structural Variants) |
| `ins_count` | Integer | Count of Insertions (Structural Variants) |
| `complex_count` | Integer | Count of Complex variants |
| `other_count` | Integer | Count of Other variant types |
| `variant_density_per_mb` | Numeric | Number of variants per Megabase |

## 2. Exclusion Metrics

**Function:** `load_exclusion_metrics()`

**Description:** Overlap metrics for exclusion regions. **Available for v5.0q benchmarks only.**

| Column | Type | Description |
| :--- | :--- | :--- |
| `bench_version` | Character | Benchmark version (e.g., `v5.0q`) |
| `ref` | Character | Reference genome |
| `bench_type` | Factor | Benchmark set type |
| `exclusions` | Character | Name of the exclusion region |
| `exclusion_bp` | Numeric | Total size of the exclusion region |
| `intersect_bp` | Numeric | Overlap with benchmark (Note: variable name in code may be `dip_intersect_bp`) |
| `pct_of_exclusion` | Numeric | Percentage of exclusion region overlapping benchmark |
| `pct_of_dip` | Numeric | Percentage of diploid genome overlapping benchmark |
| `total_variants` | Integer | Total variants in this exclusion |

*Note: Specific count columns like `snv_count`, `ins_count`, `del_count` may be present depending on variant type.*

## 3. Reference Genome Sizes

**Function:** `load_reference_sizes()`

**Description:** Information about reference genome chromosome lengths and assembly properties.

| Column | Type | Description |
| :--- | :--- | :--- |
| `ref` | Character | Reference genome name |
| `chrom` | Character | Chromosome name (standardized with "chr" prefix) |
| `length` | Integer | Total length of the chromosome |
| `ns` | Integer | Number of 'N' bases (gaps) |
| `asm_bp` | Integer | Assembled bases (`length` - `ns`) |

## 4. Variant Table

**Function:** `load_variant_table()` or `load_all_variant_tables()`

**Description:** Detailed variant-level data. These objects can be very large. Columns may vary slightly between benchmark versions, but the core columns are listed below.

| Column | Type | Description |
| :--- | :--- | :--- |
| `bench_version` | Character | Benchmark version |
| `ref` | Character | Reference genome |
| `bench_type` | Character | Benchmark set type |
| `chrom` | Character | Chromosome |
| `pos` | Integer | 1-based start position |
| `end` | Integer | End position |
| `gt` | Character | Genotype |
| `vkx` | Character | Variant class |
| `var_type` | Character | Variant type (SNV, INDEL, etc.) |
| `len_ref` | Integer | Length of reference allele |
| `len_alt` | Integer | Length of alternate allele |
| `var_size` | Integer | Size of the variant |
| `region_ids` | Character | Comma-separated IDs of regions overlapping the variant |
| `context_ids` | Character | Comma-separated IDs of genomic contexts overlapping the variant |

*Note: Additional columns like `TRF*` (Tandem Repeat Finder) or specific VCF INFO fields may be present for v5.0q.*

## 5. Genomic Context Coverage

**Function:** `load_genomic_context_coverage()`

**Description:** Base-level coverage metrics for difficult genomic contexts (e.g., how well a specific TR region is covered).

| Column | Type | Description |
| :--- | :--- | :--- |
| `context_name` | Character | Genomic context identifier |
| `chrom` | Character | Chromosome |
| `start` | Integer | Start position (0-based) |
| `end` | Integer | End position |
| `n_overlap` | Integer | Number of overlapping intervals |
| `bases_cov` | Integer | Number of bases covered |
| `ivl_len` | Integer | Length of the interval |
| `frac_cov` | Numeric | Fraction of interval covered (`bases_cov` / `ivl_len`) |

## 6. Benchmark Regions

**Function:** `load_benchmark_regions()`

**Description:**  Combined BED-like data for all benchmark regions.

| Column | Type | Description |
| :--- | :--- | :--- |
| `bench_version` | Factor | Benchmark version |
| `ref` | Factor | Reference genome |
| `bench_type` | Factor | Benchmark set type |
| `chrom` | Factor | Chromosome |
| `start` | Integer | Start position (0-based) |
| `end` | Integer | End position |
| `interval_size` | Integer | Size of the interval (`end` - `start`) |

## 7. HG002 Q100 Assembly Size

**Function:** `load_hg002q100_size()`

**Description:** Returns a single numeric value representing the total base pairs of the HG002 Q100 maternal assembly (v1.1).
