# Snakemake Workflow

This directory contains the active Snakemake pipeline for generating benchmark analysis outputs used across `analysis/`, `manuscript/`, and `docs/`.

## Structure

```text
workflow/
├── Snakefile
├── rules/
│   ├── annotation.smk
│   ├── benchmark_comparisons.smk
│   ├── chr8_synteny.smk
│   ├── common.smk
│   ├── downloads.smk
│   ├── exclusions.smk
│   ├── genomic_context_analysis.smk
│   ├── ref_processing.smk
│   └── vcf_processing.smk
├── scripts/
│   ├── annotate_old_benchmark_status.py
│   ├── combine_beds_with_id.py
│   ├── compute_coverage_table.py
│   ├── compute_exclusion_interactions.py
│   ├── count_exclusion_variants.py
│   ├── count_variants_by_genomic_context.py
│   ├── find_chr8_inversion.py
│   ├── generate_variant_parquet.py
│   ├── logging_config.py
│   └── make_chr8_figure.py
└── envs/
    ├── biotools.yaml
    ├── downloads.yaml
    ├── plotsr.yaml
    ├── python-biotools.yaml
    ├── samtools.yaml
    ├── truvari.yaml
    └── deprecated/
```

## Rule Module Summary

- `downloads.smk`: downloads benchmark VCF/BED, dip BED, references, stratifications, exclusions.
- `ref_processing.smk`: computes reference assembled base sizes.
- `vcf_processing.smk`: VCF indexing, multiallelic split, Truvari `anno svinfo`.
- `annotation.smk`: combines BED annotations and writes fully annotated VCFs.
- `genomic_context_analysis.smk`: builds genomic context coverage and variant-by-context outputs.
- `exclusions.smk`: exclusion impact tables, interaction tables, and old-only status outputs.
- `benchmark_comparisons.smk`: small variant and structural variant comparison runs.
- `chr8_synteny.smk`: optional chr8 synteny alignment/visualization workflow.
- `common.smk`: shared helper functions and target generators.

## Canonical Targets

Defined by `rule all` in `Snakefile`:

- `results/variant_tables/{benchmark}/variants.parquet`
- `results/ref_genome_sizes/{ref}_size.tsv`
- `results/genomic_context/{benchmark}/genomic_context_coverage_table.csv`
- `results/genomic_context/{benchmark}/variants_by_genomic_context.parquet`
- `results/exclusions/{benchmark}/exclusion_impact.csv`
- `results/exclusions/{benchmark}/exclusion_interactions.csv`
- `results/exclusions/{comp_id}/old_only_summary.csv`
- chr8 synteny figure targets when enabled in config

## Running

From repository root:

```bash
# Validate workflow
snakemake -n --quiet

# Execute all targets
snakemake --cores 4 --sdm conda

# Run a specific output
snakemake --cores 4 --sdm conda \
  results/variant_tables/v5.0q_GRCh38_smvar/variants.parquet
```

## Configuration

- Main config: `config/config.yaml`
- Config schema: `config/schema/config.schema.yaml`
- Test config: `config/config.test_grch38_debug.yaml`

`Snakefile` validates the loaded config against the schema at runtime.

## Notes on Environments

Environment definitions live in `workflow/envs/`; see `workflow/envs/README.md` for rationale and constraints.
