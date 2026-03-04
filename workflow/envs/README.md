# Conda Environment Layout

This directory defines the active Conda environments used by the workflow.

## Active Environments (6)

- `biotools.yaml`
  Purpose: core CLI bioinformatics tools.
  Key packages: `bcftools=1.22`, `rtg-tools=3.13`.

- `python-biotools.yaml`
  Purpose: Python processing with interval and VCF tooling.
  Key packages: `python=3.11`, `pandas`, `pyarrow`, `bcftools=1.22`, `bedtools=2.31.1`, `tabix`.

- `samtools.yaml`
  Purpose: sequence/reference utilities.
  Key packages: `samtools=1.22`, `seqkit=2.12.0`.

- `downloads.yaml`
  Purpose: remote file retrieval.
  Key package: `wget`.

- `plotsr.yaml`
  Purpose: chr8 synteny plotting workflow.
  Key packages: `minimap2`, `samtools`, `syri`, `plotsr`, `matplotlib`, `pandas<2.0`.
  Constraint: must remain isolated due to SyRI + pandas incompatibility.

- `truvari.yaml`
  Purpose: Truvari benchmark comparison and downstream table generation.
  Key packages: `Truvari==5.4.0` (pip), `bcftools=1.20`, `bedtools`, `pandas`, `pyarrow=14.0`.
  Constraint: kept separate due to `bcftools` version divergence from `biotools.yaml`.

## Deprecated Environments

Historical definitions are retained in `deprecated/`:

- `deprecated/bcftools.yaml`
- `deprecated/python.yaml`
- `deprecated/bedtools.yaml`
- `deprecated/rtg-tools.yaml`

## Validation Commands

From repository root:

```bash
# Dry-run a single env solve
mamba env create -f workflow/envs/python-biotools.yaml -n q100-env-test --dry-run

# Validate workflow wiring
make dry-run

# Execute full workflow
make run
```

## Maintenance Guidance

- Keep strict channel priorities in env files.
- Prefer pinning versions for workflow-critical tools.
- Re-run `make dry-run` after environment changes.
- If changing tool versions used in published outputs, re-run key analysis notebooks and compare outputs.
