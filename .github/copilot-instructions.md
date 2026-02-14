# Copilot Instructions for q100-variant-benchmark-paper

## 1) What this repo is

- A Snakemake-driven genomics analysis pipeline + analysis/manuscript (Quarto/R).
- Main entry: `workflow/Snakefile` (includes rules in `workflow/rules/`).
- Configuration: `config/config.yaml` validated by `config/schema/config.schema.yaml`.
- Conda envs: `workflow/envs/*.yaml` (pinned versions; channels: `conda-forge, bioconda`).

## 2) First things to look at (in order)

1. `config/config.yaml` — defines benchmark sets, references, URLs, and file hashes. Understanding names here (e.g., `v5.0q_CHM13v2.0_smvar`) is essential.
2. `workflow/Snakefile` — shows `rule all` (final outputs) and `wrapper_prefix`; also prints pipeline start/finish messages.
3. `workflow/rules/*.smk` — rule groups: `downloads.smk`, `vcf_processing.smk`, `var_tables.smk`, `var_counts.smk`, `exclusions.smk`, `strat_metrics.smk`.
4. `workflow/envs/` — check tool versions and channels used for each rule.
5. `Makefile` — quick dev tasks (`make dry-run`, `make lint`, `make test`, `make run`, `make dag`).
6. `analysis/` and `manuscript/` — Quarto/R scripts that consume `results/` CSVs for figures/tables.

## 3) Key commands (copyable)

- Setup environment:

```bash
mamba env create -f environment.yaml
mamba activate q100-smk
```

- Validate DAG (quick):

```bash
~/micromamba/envs/q100-smk/bin/snakemake -n --quiet
```

- Lint & format:

```bash
snakemake --lint
snakefmt workflow/  # or --check
```

- Run pipeline locally (conservative):

```bash
~/micromamba/envs/q100-smk/bin/snakemake --cores 4 --sdm conda
```

- Re-run a rule:

```bash
~/micromamba/envs/q100-smk/bin/snakemake --forcerun <rule_name>
```

- Generate DAG visualization:

```bash
make dag
```

## 4) Project-specific Snakemake conventions

- Every rule MUST include: `log:`, `threads:`, `resources:` (e.g., `mem_mb`), and a `conda:` directive.
- Outputs use `ensure()` wrappers to assert `non_empty=True`, `size_gt=...`, or `sha256=` for deterministic outputs.
- Download rules: use `retries: 3` and wrap outputs in `ancient()` to avoid redownloading; validate SHA256 from `config`.
- Use `wrapper_prefix` (set in `Snakefile`) when referencing Snakemake wrappers.
- For shell blocks: redirect stderr to logs (`2>> {log}` or `2>&1 | tee {log}`) and use `set -euo pipefail` patterns.

## 5) Data flow & important files

- `downloads.smk`: fetches VCF/BED and dip-beds, validates SHA256.
- `ref_processing.smk`: prepares/reference indexes (fai, dicts).
- `vcf_processing.smk`: annotates VCFs with stratifications (bcftools/bedtools steps)
- `var_tables.smk` / `var_counts.smk`: tabulate annotated variants and stratified counts.
- `results/`: final outputs consumed by `analysis/*.qmd` and `manuscript/`.

## 6) Tests & CI

- Local test target: `make test` (runs lint + format check + dry-run).
- GitHub Actions: see `.github/workflows/tests.yml` for CI steps (lint, dry-run, Python tests).

## 7) Common debugging patterns

- If outputs missing: run `snakemake --summary` and `snakemake -n --printshellcmds` to inspect what would run.
- If conda solves slowly: use `--sdm mamba` to run with mamba instead of conda.
- For runtime errors: inspect `logs/<module>/<name>.log` — rules consistently write logs.

## 8) Code style & commits

- Follow Conventional Commits: `feat:, fix:, docs:, chore:, refactor:, test:`.
- Before committing: `snakemake --lint` and `snakefmt workflow/` must pass; run `snakemake -n`.

## 9) When editing/adding rules

- Start from `workflow/rules/` templates and ensure `ensure()` validation on outputs.
- Add or update an env YAML in `workflow/envs/` (pin versions; channels: `conda-forge, bioconda`).
- Add unit-style smoke tests under `tests/fixtures` and update `tests/conftest.py` if new fixtures needed.

## 10) Short checklist for AI edits

- Read `config/config.yaml` and `workflow/Snakefile` to learn naming/wildcards.
- Add `log`, `threads`, `resources`, `conda`, and `ensure()` to new rules.
- Use `ancient()` for downloads and add SHA256 verification using `config`.
- Run `~/micromamba/envs/q100-smk/bin/snakemake -n --quiet` and `snakemake --lint` locally before proposing PRs.

---

## Analysis (Quarto/R)

- Run the Snakemake pipeline first to generate `results/` (e.g., `~/micromamba/envs/q100-smk/bin/snakemake --cores 4 --sdm conda`).
- Activate the analysis environment:
  `mamba env create -f environment.yaml` then `mamba activate q100-smk` (or use `micromamba` equivalent).
- Render Quarto documents:
  `quarto render analysis/benchmarkset_characterization.qmd` or `quarto render analysis/*.qmd`.
- Quarto caching: `analysis/*_cache/` (remove the cache folder to clear Quarto caches).
- Generated HTML and figures are saved in `analysis/<name>_cache/html/` and `analysis/<name>_files/figure-html/`.
- Analysis R packages used: `tidyverse`, `sessioninfo`, `here`, `assertthat`, `patchwork`, `gt`, `vroom`, `furrr`.
- Use `here::here()` for project-relative paths and `vroom::vroom()` for fast CSV reads.
- `analysis/*.qmd` expects `results/variant_tables/*.csv` as inputs.
