Analysis Notebook Refactoring & Pipeline Documentation Plan

 Overview

 This plan refactors three Quarto analysis notebooks to simplify data loading and improves pipeline documentation with a comprehensive output file catalog. The goal is to streamline the data loading process by using the smallest necessary output files, creating shared loading functions, and clearly documenting all pipeline outputs.

 Background

 Current Issues:

- Scattered data loading: benchmarkset_characterization.qmd loads data throughout the notebook rather than in a dedicated section
- Complex loading code: Custom parallel reading with furrr, complex tidy functions, heavy transformations during load
- No shared infrastructure: Each notebook implements its own loading logic with duplicate path parsing
- Missing documentation: No comprehensive catalog of pipeline outputs or column schemas
- Large file usage: Some analyses load full variant tables when smaller aggregated files would suffice

 Key Pipeline Outputs Identified:

 1. results/var_counts/{benchmark}/stratification_combined_metrics.csv - Primary analysis file (small, ~5-10 KB)
 2. results/exclusions/{benchmark}/exclusions_intersection_table.csv - Exclusion overlaps (v5.0q only)
 3. results/variant_tables/{benchmark}/variants.tsv - Full variant data (large, ~GB)
 4. results/diff_region_coverage/{benchmark}/{strat}_cov.bed - Base-level coverage
 5. results/ref_genome_sizes/{ref}_size.tsv - Reference genome sizes

 Benchmark Naming Pattern: {bench_version}_{ref}_{var_type} (e.g., v5.0q_GRCh38_smvar)

 ---
 Implementation Tasks

 1. Create Mermaid Diagram of Output Relationships

 File: docs/diagrams/output-relationships.mmd

 Structure:

- Show primary analysis files (small, essential) in green boxes
- Show detailed data files (large, optional) in beige boxes
- Show metadata patterns (path encoding, shared columns) in purple boxes
- Use solid arrows for direct data relationships via shared columns
- Use dashed arrows for metadata encoded in file paths

 Key relationships to visualize:

- stratification_combined_metrics.csv links to diff_region_coverage/*.bed via strat_name column
- All files link via benchmark identifier pattern in paths
- strat_name values: HP, MAP, SD, SD10kb, TR, TR10kb
- Exclusion files only exist for v5.0q benchmarks

 Groupings:

- Primary Analysis Files (stratification_combined_metrics, exclusions_intersection_table, ref_genome_sizes)
- Detailed Data Files (variant_tables, diff_region_coverage)
- Shared Metadata (path patterns, stratification names)

 ---

 1. Create Shared Data Loading Functions

 File: R/data_loading.R (new file)

 Functions to implement:

 2.1 parse_benchmark_id(file_path)

- Purpose: Extract benchmark metadata from file path
- Input: File path string or benchmark ID
- Returns: Named list with bench_version, ref, var_type
- Pattern: Extract v[0-9.]+q?_[A-Za-z0-9.]+_(smvar|stvar)
- Error handling: Stop with informative message if pattern not found

 2.2 load_stratification_metrics(results_dir, benchmark_filter)

- Purpose: Load primary analysis file with metrics + variant counts
- Input:
  - results_dir - Path to results directory (default: here::here("results"))
  - benchmark_filter - Optional vector of benchmark IDs to filter
- Returns: Tibble with columns: bench_version, ref, var_type, strat_name, strat_bp, intersect_bp, pct_of_strat, pct_of_bench, total_variants, variant counts by type, variant_density_per_mb
- Process:
   a. Find all stratification_combined_metrics.csv files
   b. Read with read_csv() using col_types "cddddiiiiiiid"
   c. Parse benchmark ID from directory name
   d. Add bench_version, ref, var_type columns
   e. Apply filter if provided

 2.3 load_exclusion_metrics(results_dir)

- Purpose: Load exclusion intersection tables (v5.0q only)
- Returns: Tibble with bench_version, ref, var_type, exclusions, exclusion_bp, intersect_bp, pct_of_exclusion, pct_of_dip
- Note: Warn if no files found (expected for non-v5.0q)

 2.4 load_reference_sizes(results_dir)

- Purpose: Load reference genome sizes
- Returns: Tibble with ref, chrom, asm_bp
- Process:
   a. Read all *_size.tsv files
   b. Calculate asm_bp = length - N
   c. Standardize chromosome names (add "chr" prefix)

 2.5 load_variant_table(benchmark_id, results_dir, filters)

- Purpose: Load full variant table (use sparingly - large file)
- Input:
  - benchmark_id - Single benchmark identifier
  - filters - List with variant_types, chromosomes, in_benchmark_only
- Returns: Tibble with variant-level data
- Warning: Display message about large file size

 2.6 load_diff_coverage(benchmark_id, results_dir, strat_filter)

- Purpose: Load base-level difficult region coverage
- Input:
  - benchmark_id - Single benchmark identifier
  - strat_filter - Optional vector of stratification names (HP, TR, etc.)
- Returns: Tibble with strat_name, chrom, start, end, n_overlap, bases_cov, ivl_len, frac_cov

 Documentation standards:

- Use roxygen2-style comments with @param, @return, @examples
- Include error handling with informative messages
- Add validation checks (file existence, expected columns)

 ---

 1. Refactor Analysis Notebooks

 3.1 benchmarkset_characterization.qmd

 Changes:

 1. Add dedicated "Data Loading" section after setup chunk (around line 270)
 2. Load primary data:

## Data Loading

# Primary metrics (small files)

 strat_metrics_df <- load_stratification_metrics()
 ref_sizes_df <- load_reference_sizes()
 exclusion_df <- load_exclusion_metrics()

# Benchmark regions

 bench_region_files <- list.files(
     path = here("resources/benchmarksets"),
     pattern = "_benchmark.bed",
     full.names = TRUE
 ) %>%
     set_names(str_extract(., "(?<=benchmarksets/).*(?=_benchmark.bed)"))

 regions_df <- bench_region_files %>%
     map_dfr(read_tsv, col_names = c("chrom", "start", "end"),
             col_types = "cii", show_col_types = FALSE, .id = "benchmark") %>%
     mutate(
         interval_size = end - start,
         chrom = if_else(str_detect(chrom, "chr"), chrom, str_c("chr", chrom)),
         bench_meta = map(benchmark, parse_benchmark_id),
         bench_version = map_chr(bench_meta, "bench_version"),
         ref = map_chr(bench_meta, "ref"),
         var_type = map_chr(bench_meta, "var_type")
     ) %>%
     select(-bench_meta, -benchmark)
 3. Add data frame documentation:

### Key Data Frames for Analysis

# 1. strat_metrics_df: Per-stratification metrics with variant counts

# Use for: variant count comparisons, density calculations

# Key columns: bench_version, ref, var_type, strat_name, variant counts

# 2. regions_df: Benchmark region intervals

# Use for: coverage calculations, interval size distributions

# Key columns: bench_version, ref, var_type, chrom, start, end, interval_size

# 3. ref_sizes_df: Reference genome sizes

# Use for: percentage calculations, normalization

# Key columns: ref, chrom, asm_bp

# 4. exclusion_df: Exclusion overlaps (v5.0q only)

# Use for: understanding what was excluded

# Key columns: bench_version, ref, var_type, exclusions, intersect_bp

 1. Remove/simplify:

- Remove tidy_smvar() and tidy_stvar() functions (lines 135-263)
- Remove complex parallel reading setup (line 267: plan(multisession, workers = 12))
- Remove read_bench_file() function
- Simplify or remove variant table loading (lines 269-298) - add comment that full variant tables should only be loaded when variant-level detail is required

 1. Update analysis sections to use strat_metrics_df instead of bench_vars_tbls where possible

 Lines to modify: 135-298 (loading section), analysis sections that reference bench_vars_tbls

 3.2 benchmark_difficult.qmd

 Changes:

 1. Update data loading section (lines 27-76):

## Data Loading

# Load difficult region coverage using shared function

 diff_cov_df <- c(
     "v5.0q_GRCh38_smvar",
     "v5.0q_GRCh37_stvar",
     "v4.2.1_GRCh38_smvar",
     "v0.6_GRCh37_stvar"
 ) %>%
     map_dfr(
         ~ load_diff_coverage(.x) %>%
             mutate(
                 bench_meta = parse_benchmark_id(.x),
                 bench_version = bench_meta$bench_version,
                 ref = bench_meta$ref,
                 var_type = bench_meta$var_type
             ) %>%
             select(-bench_meta),
         .id = "benchmark"
     )

# Calculate totals

 diff_cov_total_df <- diff_cov_df %>%
     group_by(bench_version, ref, var_type, strat_name) %>%
     summarise(
         bases_covered = sum(bases_covered),
         interval_bp = sum(interval_bp),
         .groups = "drop"
     ) %>%
     mutate(frac_cov = 1 - (bases_covered / interval_bp))

 Lines to modify: 27-76 (data loading section)

 3.3 external_evaluation.qmd

 Changes:

 1. Add section header before data loading (line 25):

## Data Loading

### External Evaluation Data

 1. No other changes - This notebook is already clean and loads external data (not pipeline outputs)

 Lines to modify: 25 (add header only)

 ---

 1. Create Pipeline Output Documentation

 4.1 Create docs/pipeline-outputs.md

 Structure:

 1. Overview - Purpose of this document, output directory structure tree
 2. File Catalog - Complete list of output types with:

- File path pattern with wildcards
- Purpose and use case
- File size estimate (small/large)
- Column schema (table format)
- Example data snippet
- Dependencies (what inputs required)
- Related outputs

 1. Column Relationships - How files link together via shared columns
 2. Metadata Encoding - How benchmark information is encoded in paths
 3. Usage Examples - Code snippets showing how to load each output type

 Primary output types to document:

 1. stratification_combined_metrics.csv - Detailed column table with 13 columns
 2. exclusions_intersection_table.csv - Column table with 5 columns
 3. variants.tsv - Column table with 10+ columns (large file warning)
 4. diff_region_coverage/*.bed - Standard bedtools coverage format (7 columns)
 5. ref_genome_sizes/*_size.tsv - Simple 4-column format

 Format for each output:

### Output Name

 **Path:** `results/category/{benchmark}/filename.ext`

 **Purpose:** Brief description of what this file contains and when to use it

 **Size:** Small (~KB) or Large (~GB)

 **Column Schema:**

 | Column | Type | Description |
 |--------|------|-------------|
 | col1 | string | Description |
 | col2 | integer | Description |

 **Usage Notes:**

- Important caveats
- Performance considerations
- When to use vs. alternative files

 **Example:**

 ```r
 # Load and use this file
 data <- load_function("benchmark_id")

 **Estimated length:** 800-1000 lines

 #### 4.2 Create `docs/data-dictionary.md`

 **Structure:**
 1. **Metrics Definitions** - Detailed explanations of each metric:
    - Definition
    - Calculation formula (mathematical notation)
    - Units
    - Interpretation guidance
    - Example values

 2. **Stratification Regions** - Description of each stratification:
    - HP, MAP, SD, SD10kb, TR, TR10kb
    - Source and version
    - Biological/technical significance
    - Typical coverage in benchmarks

 3. **Exclusion Categories** - Description of each exclusion (v5.0q only):
    - consecutive-svs, flanks, satellites, segdups, etc.
    - Definition and rationale
    - Typical exclusion size

 4. **Variant Classifications** - Definitions of variant types:
    - SNP, INDEL, DEL, INS, COMPLEX, OTHER
    - Size thresholds
    - Classification logic

 **Example format:**
 ```markdown
 ### pct_of_bench

 **Definition:** Percentage of benchmark regions that overlap with a given stratification

 **Formula:**
 pct_of_bench = (intersect_bp / total_benchmark_bp) × 100

 **Units:** Percentage (0-100%)

 **Interpretation:**
 - Higher values indicate that more of the benchmark overlaps with this difficult region
 - Values >10% suggest significant representation in the benchmark
 - Compare across benchmark versions to track coverage improvements

 **Example:**
 - HP stratification with pct_of_bench = 15.3% means that 15.3% of all benchmark bases are in homopolymer regions

 Estimated length: 600-800 lines

 4.3 Update docs/architecture.md

 Add new section after "Output Structure" (around line 80):

 ### Primary vs. Detailed Output Files

 The pipeline generates two tiers of outputs optimized for different use cases:

 **1. Primary Analysis Files** (Small, ~KB per file):
 - `stratification_combined_metrics.csv` - Aggregated metrics with variant counts
 - `exclusions_intersection_table.csv` - Exclusion overlap summaries
 - `ref_genome_sizes/*.tsv` - Reference genome metadata

 **Characteristics:**
 - Fast to load and process
 - Pre-aggregated for common analyses
 - Should be the default choice for most analyses
 - Generated by combining intermediate results

 **2. Detailed Data Files** (Large, ~MB-GB per file):
 - `variant_tables/*/variants.tsv` - Full variant-level annotations
 - `diff_region_coverage/*/*.bed` - Base-level coverage data

 **Characteristics:**
 - Slower to load and process
 - Required for variant-level or base-level investigations
 - Use with filters to reduce memory usage
 - Load only when primary files insufficient

 ### Recommended Analysis Workflow

 1. **Start with primary files** using `R/data_loading.R` functions
 2. **Generate summary statistics** and visualizations
 3. **If variant-level detail needed**, load detailed files with filters:
    ```r
    # Only load specific chromosomes or variant types
    vars <- load_variant_table(
        "v5.0q_GRCh38_smvar",
        filters = list(
            chromosomes = c("chr1", "chr2"),
            variant_types = c("SNP"),
            in_benchmark_only = TRUE
        )
    )
 4. Document data sources used for each analysis result

 **Lines to modify:** Insert new section after line 80 (Output Structure section)

 #### 4.4 Update `README.qmd`

 **Add to "Key Outputs" section** (around line 180):

 ```markdown
 ## Analysis Resources

 ### Data Loading Functions
 The `R/data_loading.R` module provides standardized functions for loading pipeline outputs:
 - `load_stratification_metrics()` - Primary analysis data with variant counts
 - `load_exclusion_metrics()` - Exclusion overlaps (v5.0q only)
 - `load_reference_sizes()` - Reference genome sizes
 - `load_variant_table()` - Full variant data (use sparingly)
 - `load_diff_coverage()` - Base-level coverage data

 See function documentation for usage examples.

 ### Output Documentation
 - **[Pipeline Outputs Reference](docs/pipeline-outputs.md)** - Complete catalog of output files with column schemas
 - **[Data Dictionary](docs/data-dictionary.md)** - Metrics definitions and calculations
 - **[Output Relationships Diagram](docs/diagrams/output-relationships.mmd)** - Visual guide to file dependencies

 Lines to modify: Add new section after line 180 (after stratification table)

 ---
 5. Verification Plan

 5.1 Function Testing

 Create test file: tests/test_data_loading.R

 Test cases:
 1. parse_benchmark_id() handles various path formats
 2. load_stratification_metrics() returns expected structure
 3. load_reference_sizes() calculates asm_bp correctly
 4. load_exclusion_metrics() handles missing files gracefully
 5. All load functions validate column schemas

 Run tests:
 library(testthat)
 source(here::here("R/data_loading.R"))
 test_file("tests/test_data_loading.R")

 5.2 Notebook Refactoring Validation

 For each refactored notebook:

 1. Add checkpoint tests in the data loading section:
 ## Data Loading Verification

 # Test 1: Expected number of benchmarks
 assertthat::assert_that(
     nrow(strat_metrics_df) > 40,  # 8 benchmarks * 6 strats
     msg = "Expected at least 48 rows in strat_metrics_df"
 )

 # Test 2: All required stratifications present
 assertthat::assert_that(
     all(c("HP", "MAP", "SD", "TR") %in% unique(strat_metrics_df$strat_name)),
     msg = "Missing expected stratification names"
 )

 # Test 3: Variant counts are non-negative
 assertthat::assert_that(
     all(strat_metrics_df$total_variants >= 0),
     msg = "Found negative variant counts"
 )
 2. Render notebook and check for errors:
 quarto render analysis/benchmarkset_characterization.qmd
 3. Compare outputs with original versions:
   - Save figures before refactoring: manuscript/figs/*_OLD.png
   - Save tables before refactoring: manuscript/tables/*_OLD.docx
   - After refactoring, visually compare or use image diff tools
   - Verify key statistics match (variant counts, region sizes, percentages)

 5.3 Documentation Validation

 1. Render mermaid diagram:
 # Using mermaid-cli
 mmdc -i docs/diagrams/output-relationships.mmd -o docs/diagrams/output-relationships.png
 2. Check markdown rendering:
   - Open docs/pipeline-outputs.md in GitHub or markdown viewer
   - Verify tables render correctly
   - Check that code blocks highlight properly
 3. Verify examples work:
   - Copy code examples from documentation
   - Run in R console with actual pipeline outputs
   - Confirm they produce expected results

 5.4 End-to-End Test

 # 1. Re-run pipeline subset to ensure outputs exist
 snakemake --cores 4 \
     results/var_counts/v5.0q_GRCh38_smvar/stratification_combined_metrics.csv \
     results/exclusions/v5.0q_GRCh38_smvar/exclusions_intersection_table.csv

 # 2. Source new data loading functions
 Rscript -e "source('R/data_loading.R'); metrics <- load_stratification_metrics(); print(dim(metrics))"

 # 3. Render all notebooks
 quarto render analysis/benchmarkset_characterization.qmd
 quarto render analysis/benchmark_difficult.qmd
 quarto render analysis/external_evaluation.qmd

 # 4. Check for rendering errors
 echo $?

 ---
 Implementation Sequence

 Recommended order to minimize risk:

 Phase 1: Foundation (1-2 days)

 1. Create R/data_loading.R with all 6 functions
 2. Create tests/test_data_loading.R
 3. Test functions thoroughly in R console
 4. Fix any issues before proceeding

 Phase 2: Documentation (1-2 days)

 1. Create docs/diagrams/output-relationships.mmd with mermaid diagram
 2. Create docs/pipeline-outputs.md with complete catalog
 3. Create docs/data-dictionary.md with metrics definitions
 4. Update docs/architecture.md (add 2 new sections)
 5. Update README.qmd (add analysis resources section)
 6. Verify all markdown renders correctly

 Phase 3: Notebook Refactoring (2-3 days)

 1. Start with benchmark_difficult.qmd (simplest, minimal changes)
   - Update data loading section
   - Add verification tests
   - Render and compare outputs
 2. Then refactor external_evaluation.qmd (add header only)
   - Minimal changes
   - Verify rendering
 3. Finally refactor benchmarkset_characterization.qmd (most complex)
   - Add data loading section
   - Remove custom functions
   - Add data frame documentation
   - Add checkpoint tests throughout
   - Render and compare all outputs (figures, tables)
   - Verify key statistics match original

 Phase 4: Verification (1 day)

 1. Run full test suite (tests/test_data_loading.R)
 2. Render all notebooks fresh
 3. Compare all figure and table outputs with originals
 4. Document any intentional differences
 5. Create verification summary report

 Phase 5: Finalization (0.5 day)

 1. Review all documentation for completeness
 2. Check for broken links in markdown
 3. Add usage examples where missing
 4. Update any outdated references
 5. Commit with clear message documenting changes

 Total Estimated Time: 5.5-8.5 days

 ---
 Critical Files Reference

 Files to create:
 - R/data_loading.R - Shared data loading functions (~350 lines)
 - docs/diagrams/output-relationships.mmd - Mermaid diagram (~100 lines)
 - docs/pipeline-outputs.md - Output file catalog (~900 lines)
 - docs/data-dictionary.md - Metrics definitions (~700 lines)
 - tests/test_data_loading.R - Unit tests (~150 lines)

 Files to modify:
 - analysis/benchmarkset_characterization.qmd (lines 135-298 + analysis sections)
 - analysis/benchmark_difficult.qmd (lines 27-76)
 - analysis/external_evaluation.qmd (line 25 - add header)
 - docs/architecture.md (add sections after line 80)
 - README.qmd (add section after line 180)

 Key pipeline outputs to document:
 - results/var_counts/{benchmark}/stratification_combined_metrics.csv
 - results/exclusions/{benchmark}/exclusions_intersection_table.csv
 - results/variant_tables/{benchmark}/variants.tsv
 - results/diff_region_coverage/{benchmark}/{strat}_cov.bed
 - results/ref_genome_sizes/{ref}_size.tsv

 Existing patterns to follow:
 - Roxygen2 documentation style from other R functions
 - ASCII diagrams in docs/architecture.md
 - Table-based documentation in docs/api-reference.md
 - Conventional commits for version control

 ---
 Success Criteria

 1. All data loading functions work correctly and have tests
 2. All three notebooks render without errors and produce equivalent outputs
 3. Mermaid diagram accurately visualizes output relationships
 4. Documentation is complete, accurate, and renders properly
 5. Primary files (stratification_combined_metrics.csv) are used instead of large variant tables where appropriate
 6. Data loading sections are clearly separated and documented in all notebooks
 7. Verification tests pass and confirm outputs match original analyses
