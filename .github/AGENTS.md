# AGENTS.md — Quick rules for AI agents

Purpose: give short, prescriptive steps so an agent can make safe, minimal-code changes.

- Look first: `config/config.yaml`, then `workflow/Snakefile`, then `workflow/rules/`.
- Always run: `~/micromamba/envs/q100-smk/bin/snakemake -n --quiet` and `snakemake --lint` locally (or ask user to run) before proposing a PR.
- For new rules:
  - Add `log`, `threads`, `resources`, and `conda` directives.
  - Wrap outputs with `ensure()`; add `sha256` for deterministic outputs when possible.
  - Add `retries: 3` and `ancient()` for download rules.
- **Terminology consistency** (when refactoring names across the pipeline):
  1. Update Snakemake rules: rule names, wildcards, `wildcard_constraints` in `common.smk`
  2. Update file paths and helper functions in `workflow/rules/common.smk`
  3. Update input/output paths in dependent rules
  4. Update Python scripts in `workflow/scripts/` that read config values or parse wildcards
  5. Update column names in all CSV/TSV outputs consistently
  6. Test with: `snakemake --lint` and dry-run on a single small benchmark
- Do not change `workflow/envs/*` without confirming tool version compatibility with project maintainers.
- **Cache management**: When modifying Snakemake rules that generate files, or updating helper scripts
  (e.g., `workflow/scripts/generate_header_lines.py`), clear affected output directories to force regeneration:
  ```bash
  rm -rf results/<target_dir>  # e.g., results/generate_annotation_headers/
  rm -rf results/<downstream_dir>  # remove dependents too
  ```
  Snakemake does not automatically detect when a script changes—cached outputs remain unchanged.
- Tests: add fixtures to `tests/fixtures` and update `tests/conftest.py` for integration-level smoke tests.
- Commit messages must use Conventional Commits. Keep changes small and focused.

Where to run tests/commands (local dev):
```bash
mamba activate q100-smk
make test
make dry-run
```

If you need to add files, place them under `workflow/` or `workflow/scripts/` and update `workflow/Snakefile` to include them.

## Analysis / Quarto / R
- Run Quarto documents after pipeline outputs exist; Quarto files live in `analysis/` and `manuscript/`.
- Use the project environment:
```bash
mamba activate q100-smk
quarto render path/to/doc.qmd
```
- Open the project with Posit/ RStudio for interactive R sessions; use the project's R package setup and avoid uncoordinated package upgrades.
- Quarto caches are safe to clear to force rebuilds; final figures and tables are expected in `analysis/` and `results/`.
- Do not commit generated HTML/PDF output unless instructed; commit source `.qmd` and code only.
