# GIAB v5q HG002 Variant Benchmark Set Analysis

## Overview

This repository contains the Snakemake workflow and analysis notebooks for the GIAB v5q HG002 variant benchmark set. The manuscript is drafted in Google Docs; this repository holds the analysis code, generated figures, and supporting documentation.

Primary goals:

1. Download and validate benchmark/reference inputs.
2. Annotate benchmark VCFs with genomic context and benchmark/exclusion regions.
3. Generate per-benchmark variant tables and genomic context summaries.
4. Quantify exclusion impact and cross-version benchmark differences.
5. Produce figures/tables used in the manuscript.

## Current Workspace Layout

```text
project-root/
├── workflow/                 # Snakemake pipeline (rules, scripts, envs)
├── config/                   # Runtime config + JSON schema
├── R/                        # Shared R helpers (loading, cache, themes)
├── analysis/                 # Quarto analysis notebooks
├── manuscript/figs/           # Generated figures for manuscript
├── docs/                     # Architecture, outputs, troubleshooting
├── notes/                    # Dated working notes, plans, and prep checklists
├── tests/                    # R + Python tests
├── scripts/                  # Utility scripts
├── resources/                # Downloaded pipeline inputs (gitignored)
├── results/                  # Pipeline outputs (gitignored)
└── logs/                     # Workflow logs (gitignored)
```

## Pipeline Entry Points

- Workflow entry: `workflow/Snakefile`
- Main workflow config: `config/config.yaml`
- Config schema: `config/schema/config.schema.yaml`
- Developer commands: `Makefile`

The active rule modules included by `workflow/Snakefile` are:

- `rules/common.smk`
- `rules/downloads.smk`
- `rules/ref_processing.smk`
- `rules/vcf_processing.smk`
- `rules/annotation.smk`
- `rules/genomic_context_analysis.smk`
- `rules/exclusions.smk`
- `rules/benchmark_comparisons.smk`
- `rules/chr8_synteny.smk`

## Key Outputs

- `results/variant_tables/{benchmark}/variants.parquet`
- `results/genomic_context/{benchmark}/genomic_context_coverage_table.csv`
- `results/genomic_context/{benchmark}/variants_by_genomic_context.parquet`
- `results/exclusions/{benchmark}/exclusion_impact.csv`
- `results/exclusions/{benchmark}/exclusion_interactions.csv`
- `results/exclusions/{comp_id}/old_only_summary.csv`
- `results/ref_genome_sizes/{ref}_size.tsv`
- `results/chr8_synteny/chr8_figure.{pdf,png}` (when chr8 synteny is enabled)

## Getting Started

## Prerequisites

- Conda or Mamba
- Snakemake >= 8
- Quarto (for notebook/manuscript/site rendering)
- R with project packages (renv recommended)

## Environment setup

```bash
mamba env create -f environment.yaml
mamba activate q100-smk
```

## Run workflow

```bash
# Validate DAG and inputs
make dry-run

# Run workflow
make run

# Lint + format checks + dry-run validation
make test
```

## Render notebooks/manuscript

```bash
# Example notebook
quarto render analysis/benchmarkset_characterization.qmd

# Project root Quarto document
quarto render index.qmd
```

## Build site

The repository is configured as a Quarto website project. To render:

```bash
quarto render           # Render full site to _site/
quarto preview          # Live preview with auto-reload
quarto publish gh-pages # Publish to GitHub Pages
```

## Documentation Index

- `docs/architecture.md`: system architecture and design rationale
- `docs/pipeline-outputs.md`: output file catalog
- `docs/data-dictionary.md`: metric/column definitions
- `docs/api-reference.md`: shared helper/API docs
- `docs/troubleshooting.md`: common failures and fixes
- `notes/README.md`: index of dated/working notes
- `analysis/README.md`: notebook-specific guidance
- `workflow/README.md`: workflow module and script guide
- `tests/README.md`: test execution and fixture details

## Artifact Policy

- Commit source and durable assets: workflow code (`workflow/`, `config/`, `R/`, `scripts/`), notebooks (`analysis/*.qmd`), and docs.
- Do not commit generated run outputs: `results/`, `resources/`, `logs/`, `.snakemake/`, `analysis/cache/`, rendered notebook HTML, and site build outputs under `_site/` and `_freeze/`.
- Keep root-level exports/archive files out of git: LaTeX export intermediates (`index.tex`), run reports (`pipeline_run.html`), and archive bundles (`*.tar.gz`).
- Generated figures for the manuscript are versioned under `manuscript/figs/`.

## Workspace Notes

- Large/generated directories (`results/`, `resources/`, `logs/`, `.snakemake/`, `analysis/cache/`, rendered HTML) are expected in local runs and are mostly gitignored.
- Manuscript figures under `manuscript/figs/` are versioned.

## License

This project is licensed under CC0 1.0 Universal.
