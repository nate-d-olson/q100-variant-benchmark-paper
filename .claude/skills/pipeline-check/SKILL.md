---
name: pipeline-check
description: Run Snakemake linting, dry-run validation, and R linting to validate pipeline changes
disable-model-invocation: true
---

# Pipeline Check

Run the full pipeline validation suite to catch issues before committing.

## Steps

1. **Snakemake lint**: Run `micromamba run -n q100-smk snakemake --lint` to check for rule issues
2. **Snakemake dry-run**: Run `micromamba run -n q100-smk snakemake -n -- all` to verify the DAG resolves
3. **Snakefmt check**: Run `micromamba run -n q100-smk snakefmt --check workflow/` to verify Snakemake formatting
4. **R lint**: Run `Rscript -e 'lintr::lint_dir("R/")'` to check R code quality
5. **Markdownlint**: Run `markdownlint '**/*.md' --ignore node_modules` to check markdown files

Report results clearly, grouping errors by category (lint, dry-run, format, R, markdown). If everything passes, confirm all checks passed.
