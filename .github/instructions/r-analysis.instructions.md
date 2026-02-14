---
applyTo: "R/**/*.R,analysis/**/*.{R,qmd}"
---

# R and Quarto Analysis Instructions

This document provides guidelines for working with R code and Quarto analysis documents in this project.

## Code Style and Formatting

**Formatter:** `air` (Automatic Indentation for R)
- Configuration: `air.toml` in project root
- Line width: 100 characters
- Indent: 2 spaces
- Format command: `air format .`
- Check without modifying: `air format --check .`

**Linter:** `lintr`
- Configuration: `.lintr` in project root
- Allows `dotted.case` for internal helper functions
- Run linting: `Rscript -e "lintr::lint_dir()"`

**Pipe Operator:** Use `%>%` (magrittr pipe) instead of `|>` (native pipe)
- Reason: Consistency with existing codebase
- The `data_loading.R` defines `%>%` if missing for standalone usage

## Data Loading Infrastructure

The project uses a sophisticated caching system for data loading:

### Key Files
- **R/data_loading.R** - Main loading functions (sources schemas.R and cache.R)
- **R/schemas.R** - Arrow schema registry, factor levels, validation rules
- **R/cache.R** - Parquet caching with zstd compression
- **docs/data-dictionary.md** - Detailed documentation of all data objects

### Caching System
- **Format:** Parquet with zstd compression (level 3) via Arrow
- **Location:** `analysis/cache/` (configurable via `getOption("q100.cache_dir")`)
- **Invalidation:** Cache key = `rlang::hash()` of dataset name + source file mtimes + params
- **Metadata:** Pipeline metadata stored as JSON in Parquet file-level key-value metadata
- **Factors:** Stripped to character before write, restored on read from schema registry

### Loading Functions

All loading functions follow this pattern:
```r
load_<dataset> <- function(use_cache = TRUE, force_refresh = FALSE, ...) {
  # Check cache first
  if (use_cache && !force_refresh) {
    cached <- read_cache("dataset_name", ...)
    if (!is.null(cached)) return(cached)
  }
  
  # Load data
  data <- # ... load logic ...
  
  # Write cache (wrapped in tryCatch for safety)
  if (use_cache) {
    write_cache(data, "dataset_name", ...)
  }
  
  return(data)
}
```

### Key Loading Functions
- `load_genomic_context_metrics()` - Primary analysis data with variant counts per genomic context
- `load_primary_analysis_data()` - Load commonly used validated data frames in one call
- `load_exclusion_metrics()` - Exclusion overlaps (v5.0q only)
- `load_variant_table()` - Full variant data (large files, use sparingly)
- `load_genomic_context_coverage()` - Base-level coverage data
- `parse_benchmark_id()` - Parse benchmark metadata from file paths

## Factor Levels

Factor levels are centralized in `R/schemas.R`:
- `BENCH_VERSION_LEVELS` - Benchmark version ordering
- `CHROM_LEVELS` - Chromosome ordering
- `VAR_TYPE_LEVELS` - Variant type ordering
- `SZBIN_LEVELS` - Size bin ordering

Always use these constants when creating factors to ensure consistency.

## File Path Handling

Use `here::here()` for project-relative paths:
```r
# Good
data_path <- here::here("results", "variant_tables", "v5.0q_GRCh38_smvar.csv")

# Bad (don't use relative paths from file location)
data_path <- "../../results/variant_tables/v5.0q_GRCh38_smvar.csv"
```

## Fast CSV Reading

Use `vroom::vroom()` for fast CSV/TSV reading:
```r
# Fast reading with automatic type detection
data <- vroom::vroom(file_path)

# With column types
data <- vroom::vroom(file_path, col_types = cols(
  chrom = col_character(),
  pos = col_integer(),
  var_type = col_factor(levels = VAR_TYPE_LEVELS)
))
```

## Quarto Documents

### Setup Block Pattern

All analysis notebooks should source the setup file:
```r
#| label: setup
#| include: false

source(here::here("analysis", "_notebook_setup.R"))
source(here::here("R", "data_loading.R"))
```

### Caching

Quarto notebooks have their own caching (separate from R cache):
- Cache location: `analysis/<notebook>_cache/`
- Clear cache: `rm -rf analysis/<notebook>_cache/`
- Rendered output: `analysis/<notebook>_files/figure-html/`

### Rendering

```bash
# Single notebook
quarto render analysis/benchmarkset_characterization.qmd

# All notebooks
quarto render analysis/*.qmd

# Preview during development
quarto preview analysis/benchmarkset_characterization.qmd
```

## Testing

R tests use `testthat`:
```bash
# Run all tests
Rscript -e 'testthat::test_dir("tests")'

# Run specific test file
Rscript -e 'testthat::test_file("tests/test_cache.R")'
```

## Common Patterns

### Adding a New Cached Dataset

1. Add schema to `R/schemas.R`:
   - `get_arrow_schema()` - Arrow schema definition
   - `get_factor_levels()` - Factor level mappings
   - `get_validation_rules()` - Validation assertions

2. Add loading function to `R/data_loading.R`:
   - Add `use_cache` and `force_refresh` parameters
   - Call `read_cache()` before loading
   - Call `write_cache()` after loading
   - Wrap cache writes in `tryCatch` for safety

3. Document in `docs/data-dictionary.md`:
   - Add entry with column descriptions
   - Include example usage
   - Note any special considerations

### File Matching Pattern

When using `fs::dir_ls()`, remember that `glob` matches against **full path**:
```r
# Use grepl() on fs::path_file() for filename-only matching
files <- fs::dir_ls("results/variant_tables/", regexp = "\\.parquet$")
files <- files[grepl("v5\\.0q", fs::path_file(files))]
```

## Dependencies

Key R packages used:
- **tidyverse** - Data manipulation and visualization
- **here** - Project-relative paths
- **vroom** - Fast CSV reading
- **arrow** - Parquet I/O
- **assertthat** - Validation
- **patchwork** - Plot composition
- **gt** - Table formatting
- **furrr** - Parallel processing
- **sessioninfo** - Session information

## Package Management

- This is NOT a formal R package (no DESCRIPTION/NAMESPACE)
- `@export` roxygen tags are for documentation only
- Source files directly with `source()`
- Do not use `library()` or `require()` on project code

## Before Committing

1. Format R code: `air format .`
2. Check formatting: `air format --check .`
3. Run tests: `Rscript -e 'testthat::test_dir("tests")'`
4. Verify notebooks render: `quarto render <file>.qmd`

## Known Issues

- `parse_benchmark_id()` regex captures `"5.0q"` not `"v5.0q"` (tests expect "v" prefix)
- Some `test_data_loading.R` tests reference old column names from pre-refactor
- Tests requiring real data fail without `results/` directory (pipeline outputs not in repo)

## Maintenance Note

⚠️ **The data dictionary document (`docs/data-dictionary.md`) MUST be updated whenever changes are made to:**
- `R/data_loading.R` - Loading functions or return values
- `R/schemas.R` - Schema definitions or validation rules  
- Pipeline output formats - Column names or data types
