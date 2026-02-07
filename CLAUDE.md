# Q100 Variant Benchmark Paper

## Project Overview

Quarto manuscript analyzing the GIAB Q100 HG002 variant benchmark. The Snakemake pipeline processes variant calls across multiple reference genomes (GRCh37, GRCh38, CHM13v2.0) and benchmark versions (v0.6, v4.2.1, v5.0q).

## Key Directories

- `R/` - R source files (data loading, schemas, caching)
- `analysis/` - Quarto notebooks and cached data
- `config/` - Pipeline configuration (config.yaml)
- `docs/` - Architecture docs, data dictionary, troubleshooting
- `tests/` - R tests (testthat) and Python tests (pytest)
- `workflow/` - Snakemake rules and Python scripts
- `results/` - Pipeline outputs (gitignored, must be generated)
- `resources/` - Benchmark BED files (gitignored, must be generated)

## R Data Loading Infrastructure

### File Organization

```
R/
├── schemas.R       # Arrow schema registry, factor levels, validation rules
├── cache.R         # Parquet caching (sources schemas.R)
└── data_loading.R  # All loading functions (sources schemas.R and cache.R)
```

Source order: `data_loading.R` sources `schemas.R` and `cache.R` at the top. Use `source(here::here("R/data_loading.R"))` to get everything.

### Caching System

- **Format:** Parquet with zstd compression (level 3) via Arrow
- **Location:** `analysis/cache/` (configurable via `getOption("q100.cache_dir")`)
- **Invalidation:** Cache key = `rlang::hash()` of dataset name + source file mtimes + params
- **Metadata:** Pipeline metadata (R version, package versions, config) stored as JSON in Parquet file-level key-value metadata under key `q100_pipeline_metadata`
- **Factors:** Stripped to character before Parquet write, restored on read from schema registry

### Adding a New Cached Dataset

1. `R/schemas.R`: Add entries to `get_arrow_schema()`, `get_factor_levels()`, `get_validation_rules()`
2. `R/data_loading.R`: Add `use_cache`/`force_refresh` params, call `read_cache()` before loading, `write_cache()` after

### Testing

```bash
Rscript -e 'testthat::test_file("tests/test_cache.R")'      # Schema + cache tests (45 tests)
Rscript -e 'testthat::test_file("tests/test_data_loading.R")'  # Data loading tests (has pre-existing failures)
```

Cache tests use `withr::local_options(q100.cache_dir = tempdir)` for isolation.

### Known Issues

- `parse_benchmark_id()` regex captures `"5.0q"` not `"v5.0q"` -- tests expect `"v"` prefix
- `test_data_loading.R` references old column names (`strat_name`, `strat_filter`) from pre-refactor
- Tests that load real data fail without `results/` directory (pipeline outputs not in repo)

## Code Quality

- **Formatter:** `air` (config in `air.toml`: 100-char lines, 2-space indent)
- **Linter:** `lintr` (config in `.lintr`: allows `dotted.case` for internal helpers)

## Important Patterns

- `fs::dir_ls` glob matches against **full path** via `glob2rx` -- use `grepl()` on `fs::path_file()` for filename-only matching
- Factor levels are centralized in `R/schemas.R` constants (e.g., `BENCH_VERSION_LEVELS`, `CHROM_LEVELS`)
- All loading functions have `use_cache = TRUE` and `force_refresh = FALSE` defaults; cache writes are wrapped in `tryCatch` so failures don't break loading
