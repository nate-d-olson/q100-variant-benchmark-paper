# Q100 Variant Benchmark — Project-Specific Guidance

## Pipeline Organization

**Terminology:**
- **Genomic Context**: Difficult regions (HP, TR, SD, MAP, etc.)—used for variant analysis and metrics computation
- **Stratification**: Same regions as genomic context, but when used in GIAB comparison analysis or benchmarkset comparisons
- **Benchmark Regions**: The set of regions where variants are confident—separate from genomic contexts
- **Bench Type**: Classification for benchmark variants (smvar = small variants, stvar = structural variants)

**Key Rules by Layer:**
1. **Metrics** (`workflow/rules/strat_metrics.smk`): Compute overlap statistics between genomic contexts and benchmarks
   - Outputs: `results/genomic_context_metrics/{benchmark}/genomic_context_coverage_table.csv`
2. **Variant Tables** (`workflow/rules/var_tables.smk`): Annotate VCFs with genomic context IDs and generate TSV
   - Key rule: `annotate_vcf_genomic_contexts` → adds INFO/CONTEXT_IDS to VCF
   - Outputs: `results/variant_tables/{benchmark}/variants.tsv`
3. **Variant Counts** (`workflow/rules/var_counts.smk`): Count variants per genomic context
   - Outputs: `results/var_counts/{benchmark}/genomic_context_combined_metrics.csv`
4. **Comparisons** (`workflow/rules/comparisons.smk`): Truvari comparison between benchmark versions
   - Uses "stratification" terminology (not genomic_context) for GIAB comparison analysis

**Output File Patterns:**
- Genomic context files: `results/genomic_context_*/{benchmark}/*.csv`
- Variant annotations: `results/variant_tables/{benchmark}/*.tsv`
- Exclusions (v5.0q only): `results/exclusions/{benchmark}/*.csv`

## Common Debugging Patterns

**bcftools annotation error**: "The INFO tag 'CONTEXT_IDS' is not defined"
- **Cause**: Header file was cached with old definitions before script changes
- **Root cause**: `generate_header_lines.py` was updated but `results/generate_annotation_headers/` had stale cache
- **Fix**: Delete affected output directories to force regeneration
  ```bash
  rm -rf results/generate_annotation_headers/
  rm -rf results/annotate_vcf_genomic_contexts/
  rm -rf results/annotate_vcf_regions/
  rm -rf results/combine_genomic_context_beds/
  rm -rf results/extract_info_fields/
  rm -rf results/variant_tables/
  ```
- **Prevention**: Always clear downstream outputs when modifying scripts that generate headers/annotations

**Wildcard constraint errors**: Check `workflow/rules/common.smk` for `wildcard_constraints`
- Must match all wildcards used in affected rules
- Format: `strat_name="|".join(sorted(_all_strat_names))`
- For new genomic context rules: use `get_genomic_context_ids()` helper function

**Script changes not taking effect**:
- Snakemake does not automatically detect when a script changes
- Solution: Delete the output directory that depends on the script; re-run to regenerate
- Example: modify `combine_beds_with_id.py` → delete `results/combine_genomic_context_beds/`

## Data Loading (R/Quarto)

**Key functions in `R/data_loading.R`:**
- `load_stratification_metrics(results_dir, benchmark_filter)` — Primary analysis file
  - Returns: metrics + variant counts per genomic context
  - Use this for most analyses instead of full variant tables
  - Outputs: `results/genomic_context_metrics/{benchmark}/genomic_context_combined_metrics.csv`
- `load_exclusion_metrics(results_dir)` — Exclusion overlap tables (v5.0q only)
  - Warns if no files found
- `load_variant_table(benchmark_id, results_dir, filters)` — Full variant data
  - Large file (~GB); use filters to reduce memory
  - Only load when variant-level detail is required
- `load_diff_coverage(benchmark_id, results_dir, strat_filter)` — Base-level coverage
  - From: `results/genomic_context_coverage/{benchmark_version}_{ref}_{bench_type}/{genomic_context}_cov.bed`

**Quarto notebooks:**
- `analysis/benchmarkset_characterization.qmd` — Primary analysis; loads `stratification_metrics_df`
- `analysis/benchmark_difficult.qmd` — Coverage analysis; loads `diff_cov_df`
- `analysis/external_evaluation.qmd` — External benchmark comparisons

## Column Naming Conventions

**Genomic Context Files** (CSV/TSV):
- `context_name` — Genomic context identifier (HP, TR, SD, MAP, etc.)
- `context_bp` — Total size of genomic context
- `intersect_bp` — Overlap with benchmark
- `pct_of_context` — % of genomic context covered by benchmark
- `pct_of_bench` — % of benchmark within this genomic context
- `total_variants` — Total variant count
- `variant_density_per_mb` — Variants per Mb of benchmark-covered region

**VCF Annotation Fields** (after bcftools annotate):
- INFO/CONTEXT_IDS — Comma-separated genomic context region IDs overlapping variant
- INFO/REGION_IDS — Comma-separated benchmark + exclusion region IDs

**Variant Table** (`variants.tsv`):
- Column 9: %INFO/CONTEXT_IDS — Genomic contexts (extracted during bcftools query)
- Column 10: %INFO/REGION_IDS — Benchmark + exclusions (extracted during bcftools query)

## Benchmark ID Format

Pattern: `{bench_version}_{ref}_{bench_type}`

Examples:
- `v5.0q_GRCh38_smvar` — v5.0q small variants on GRCh38
- `v5.0q_GRCh37_stvar` — v5.0q structural variants on GRCh37
- `v4.2.1_GRCh38_smvar` — v4.2.1 small variants on GRCh38

File paths encode this pattern:
```
results/genomic_context_metrics/v5.0q_GRCh38_smvar/
results/variant_tables/v5.0q_GRCh38_smvar/variants.tsv
results/exclusions/v5.0q_GRCh38_smvar/  # only if exclusions configured
```
