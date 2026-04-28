# Troubleshooting Guide

Common pipeline and analysis issues encountered in this repo, with verified fixes.
Focus: things that are non-obvious from the error message or that have bitten us
more than once. Generic Snakemake/conda issues are not duplicated here.

## Conda / Environment Setup

### Anaconda CDN blocked by org firewall

**Symptom:** `mamba` / `conda` fails when downloading from `conda.anaconda.org` or
`repo.anaconda.com`.

**Cause:** Org firewall blocks Anaconda CDN URLs.

**Fix:** Configure `~/.condarc` to remap `conda-forge` and `bioconda` to prefix.dev:

```yaml
custom_multichannels:
  conda-forge:
    - https://prefix.dev/conda-forge
  bioconda:
    - https://prefix.dev/bioconda
default_channels: []
```

`channel_alias` does **not** work for prefix.dev — its URL path (`/conda-forge`)
differs from the conda channel-alias convention (`/channels/conda-forge`). Use
`custom_multichannels`.

Verify:

```bash
mamba search --channel bioconda samtools  # should show "prefix.dev"
```

### Snakemake-managed conda envs out of date

After editing `workflow/envs/*.yaml`, force a rebuild:

```bash
rm -rf .snakemake/conda
snakemake --cores 4 --sdm conda --conda-create-envs-only
```

The pipeline uses `--sdm conda` (Snakemake 8 syntax), not the older `--use-conda`.
See `Makefile` `make run` for the canonical invocation.

## Pipeline Execution

### Stale annotation header — "INFO tag CONTEXT_IDS not defined"

**Cause:** Header file in `results/generate_annotation_headers/` was cached before
a change to the heredoc in `workflow/rules/annotation.smk`.

**Fix:** Delete the stale outputs so they rebuild:

```bash
rm -rf results/generate_annotation_headers/
rm -rf results/annotate_vcf_genomic_contexts/
rm -rf results/annotate_vcf_regions/
rm -rf results/combine_genomic_context_beds/
rm -rf results/variant_tables/
```

Snakemake does not detect script-only changes — when modifying a script that
generates headers or annotations, always clear downstream outputs manually.

### Wildcard constraint errors

Constraints live in `workflow/rules/common.smk` under `wildcard_constraints`.
The pattern is `strat_name="|".join(sorted(_all_strat_names))`. New context rules
should use `get_stratifications_for_ref()`. Internal helpers use a `_` prefix
(e.g., `_get_exclusion_config`).

### Missing chromosome in output

Chromosome lists live in `config.yaml` under `chromosomes:`, not hardcoded.
`get_chromosomes()` reads them and applies the GRCh37/GRCh38 prefix. If a chromosome
is missing from outputs, check the config.

## Variant Tables

### Empty variant counts with malformed context IDs

**Symptom:** Logs show `Contexts: ["'MAP'", "'MAP')", "('HP'", ...]` — i.e.,
fragments of Python tuple `repr()` instead of clean context names.

**Status:** Fixed February 2026 in commit `78747aa` ("convert Truvari tuple
objects to comma-separated strings in variant tables"). The fix lives in
`normalize_annotation()` at `workflow/scripts/generate_variant_parquet.py:66`.

**Verification:** Logs should show clean names:

```
Contexts: ['HP', 'MAP', 'SD', ...]
```

If you regenerate variant Parquet from an old branch and see the broken output,
update to current `main` and re-run; no manual workaround is needed on current code.

### Comma-separated context_ids not exploded for filtering

**Symptom:** Filtering `variants_df` by a single context misses variants whose
`context_ids` contains multiple comma-separated entries.

**Status:** Fixed in commit `35bc5fd` for the fold-change analysis. When writing
new analysis code, `tidyr::separate_rows(context_ids, sep = ",")` (or the parquet
loader's filter API) before grouping by context.

## Chr8 Synteny

### SyRI crash — `ValueError: buffer source array is read-only`

**Cause:** pandas 2.0 Copy-on-Write makes DataFrame slice arrays non-writeable;
SyRI's Cython code (`synsearchFunctions.pyx:534`) needs a writable memoryview.

**Fix:** `workflow/envs/plotsr.yaml` pins `pandas<2.0`. If you see this error,
rebuild the env and clear failed outputs:

```bash
rm -rf results/chr8_synteny/syri/
snakemake chr8_synteny --sdm conda
```

### SyRI `--prefix` deprecation

**Symptom:** "For specifying output folder use --dir, use --prefix for modifying
the output file names" — may crash on some SyRI versions.

**Status:** `chr8_syri` rule already uses `--dir {params.outdir} --prefix {params.prefix}`.

### plotsr produces wrong figure with 3 genomes

**Cause:** plotsr requires *consecutive-genome* SyRI files. For layout
[REF, MAT, PAT] it needs `ref_matsyri.out` and **`mat_patsyri.out`** — not
`ref_patsyri.out`.

**Status:** Pipeline produces all three (`ref_mat`, `mat_pat`, `ref_pat`).
`chr8_make_figure` consumes `mat_patsyri.out`; `chr8_find_inversion` separately
consumes `ref_patsyri.out` to detect PAT inversions in REF coordinate space.

### `find_chr8_inversion.py` — "ValueError: int('chr8')"

**Cause:** SyRI `.out` query columns are 6 and 7 (0-indexed), not 5 and 6.
Column 5 is `qryChr` (a string).

**Status:** Fixed in `find_chr8_inversion.py`. If you write new SyRI parsers:
0=refChr, 1=refStart, 2=refEnd, 5=qryChr, 6=qryStart, 7=qryEnd, 10=type.

## Tooling

### `gh` CLI fails with x509 certificate error

**Symptom:** `gh` returns `tls: failed to verify certificate: x509: OSStatus -26276`
intermittently.

**Cause:** Org network proxy.

**Workaround:** Use the MCP GitHub tools (e.g., `mcp__plugin_github_github__merge_pull_request`)
for PR operations. Git over SSH is unaffected.

## R Caching

### Stale cache after pipeline re-run

The cache invalidates on source-file mtime changes. If a notebook still shows old
data:

```r
source(here::here("R/data_loading.R"))

variants <- load_variant_table("v5.0q_GRCh38_smvar", force_refresh = TRUE)
# or:
clear_cache()
```

### Cache write failures (non-fatal)

Loaders return data even when caching fails (warning emitted). Common causes:
disk full, write permissions on `analysis/cache/`, or schema validation failure
against `R/schemas.R`. Inspect with `cache_info()`.

### Custom cache directory

```r
options(q100.cache_dir = "/path/to/cache")
source(here::here("R/data_loading.R"))
```

## Diagnostic Commands

```bash
# Show pipeline status
snakemake --summary

# Why a rule needs to re-run
snakemake -n --reason {rule_name}

# Print shell commands for a target
snakemake --dry-run --printshellcmds {target}

# DAG / rule graph
make dag                                 # results/dag/pipeline.pdf
snakemake --rulegraph | dot -Tpng > rulegraph.png

# Lint everything (Snakemake + Python + R)
make lint

# Pre-merge gate
make test                                # lint + format-check + dry-run
```

## Quick Reference

| Symptom | Likely cause | Fix |
|---|---|---|
| "INFO tag CONTEXT_IDS not defined" | Stale annotation header cache | `rm -rf results/generate_annotation_headers/` and downstream |
| Mojibake context names in logs | Pre-`78747aa` Truvari output | Update to current main |
| SyRI `read-only buffer` crash | pandas ≥ 2.0 in plotsr env | Rebuild env (already pinned) |
| Empty `load_exclusion_metrics()` | Not a v5.0q benchmark | Expected; warning is informational |
| `gh` x509 error | Org proxy | Use MCP GitHub tools |
| Wildcard constraint mismatch | Missing entry in `common.smk` | Add to `wildcard_constraints` |
| Snakemake doesn't re-run after script edit | Snakemake doesn't track scripts | `rm -rf` the affected output dir |

## Related Docs

- [Architecture Overview](architecture.md)
- [Pipeline Outputs Reference](pipeline-outputs.md)
- `workflow/envs/README.md` — env consolidation notes
- `CLAUDE.md` — project-specific conventions
