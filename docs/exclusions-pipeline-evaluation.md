# Exclusions vs. Genomic Context Pipelines

**Status:** Decision record — explains why the exclusions pipeline (v5.0q only)
is structured differently from the genomic context pipeline. The two-pipeline
split is intentional and is not slated to converge.

> **Historical note:** An earlier version of this document (Feb 2026) evaluated
> whether `compute_bed_metrics.py` should be merged with `compute_coverage_table.py`.
> That script has since been removed entirely — `compute_exclusion_metrics` now
> uses inline `bedtools` shell commands. The architectural reasoning below still
> applies.

## Pipeline Shapes

### Genomic Context

```
resources/stratifications/{ref}_{context}.bed.gz
  ↓ genomic_context_coverage          (bedtools coverage)
results/genomic_context/{benchmark}/coverage/{context}_cov.bed
  ↓ compute_genomic_context_coverage_table   (compute_coverage_table.py)
results/genomic_context/{benchmark}/genomic_context_coverage_table.csv
```

Two-step: per-context bedtools coverage → single Python aggregator that
processes ALL contexts in one pass.

### Exclusions (v5.0q only)

```
resources/exclusions/{benchmark}/{exclusion}_{idx}.bed
  ↓ materialize_exclusion             (bedtools sort/merge — handles single + pair types)
results/exclusions/{benchmark}/{exclusion}.bed
  ↓ compute_exclusion_metrics         (inline bedtools shell)
results/exclusions/{benchmark}/coverage/{exclusion}.tsv
  ↓ compute_exclusion_impact          (count_exclusion_variants.py)
results/exclusions/{benchmark}/exclusion_impact.csv
```

Three-step: materialize → BED metrics → join with variant counts.

## Why They're Different

| Aspect | Genomic Context | Exclusions | Why |
|---|---|---|---|
| Pre-computed coverage | Yes (`*_cov.bed`) | No (direct from BED) | Coverage BEDs are reused by other notebook analyses |
| Variant integration | No (separate pipeline) | Yes (in `compute_exclusion_impact`) | Exclusion analysis specifically asks "which variants are removed?" |
| Aggregation | Single rule, all contexts at once | Per-exclusion TSV → join | Variant-level join is per-exclusion by construction |
| Input prep | Direct from `stratifications/` | `materialize_exclusion` (sort/merge) | Exclusions arrive as multiple paired files for some types |
| Conda env | `python-biotools` (single) | `python-biotools` (BED) + `truvari` (variant counts) | Variant counting uses Truvari |

`materialize_exclusion` is **not** a symlink rule — it sorts and (for `pair`-type
exclusions) merges multiple source BEDs. Treating it as a no-op would break the
downstream metrics.

## Cross-version Outputs

The exclusions pipeline also produces `results/exclusions/{comp_id}/old_only_*.csv`
files via `annotate_old_benchmark_status` for cross-version analysis (Q3:
"what regions/variants were in v4.2.1 or v0.6 but excluded from v5.0q?").

## When Would We Reconsider?

A merge would be worth revisiting if:

- Exclusion BED metrics become the runtime bottleneck (currently dominated by
  variant counting, which is well-isolated in `count_exclusion_variants.py`)
- A new analysis needs per-exclusion coverage BEDs the way genomic-context
  analyses use `*_cov.bed`
- The exclusion set grows large enough that the per-exclusion TSV materialization
  becomes a maintenance burden

None of these apply today. The pipelines are appropriately different — not
inconsistent.

## Related

- [Architecture Overview](architecture.md)
- [Pipeline Outputs Reference](pipeline-outputs.md)
- `workflow/rules/exclusions.smk` — current implementation
- `workflow/rules/genomic_context_analysis.smk` — counterpart
