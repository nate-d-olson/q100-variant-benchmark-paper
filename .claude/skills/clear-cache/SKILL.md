---
name: clear-cache
description: Clear stale pipeline caches based on which script or rule file changed
disable-model-invocation: true
---

# Clear Cache

When a pipeline script or rule file changes, Snakemake may not detect the change automatically. This skill identifies which output directories depend on the changed file and removes them to force regeneration.

## Usage

Invoke with an optional argument specifying the changed file path. If no argument is given, check `git diff --name-only` for recently modified workflow scripts.

## Dependency Map

Use these known dependency chains to determine what to clear:

### Header/Annotation scripts

If `generate_header_lines.py`, `combine_beds_with_id.py`, or annotation-related scripts changed:

```
rm -rf results/generate_annotation_headers/
rm -rf results/annotate_vcf_genomic_contexts/
rm -rf results/annotate_vcf_regions/
rm -rf results/combine_genomic_context_beds/
rm -rf results/variant_tables/
```

### Variant table scripts

If `generate_variant_parquet.py` or variant classification scripts changed:

```
rm -rf results/variant_tables/
rm -rf results/var_counts/
```

### Variant count scripts

If `count_variants_by_genomic_context.py` changed:

```
rm -rf results/var_counts/
```

### Metrics scripts

If `compute_bed_metrics.py` or metrics scripts changed:

```
rm -rf results/genomic_context_metrics/
```

### Exclusion scripts

If exclusion-related scripts changed:

```
rm -rf results/exclusions/
```

### R cache

If `R/schemas.R`, `R/cache.R`, or `R/data_loading.R` changed:

```
rm -rf analysis/cache/
```

## Rules

1. Always show the user which directories will be deleted and get confirmation before proceeding
2. After clearing, suggest running `micromamba run -n q100-smk snakemake -n -- all` to verify the DAG
3. If the changed file doesn't match any known dependency, warn the user and suggest checking the Snakemake rule that uses that script
