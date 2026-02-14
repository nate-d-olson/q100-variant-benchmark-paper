---
applyTo: "**"
---

# Development Workflow Guide

This document describes the development workflow for the GIAB v5q benchmark set manuscript.

## Environment Setup

### Initial Setup

```bash
# Clone repository
git clone https://github.com/nate-d-olson/q100-variant-benchmark-paper.git
cd q100-variant-benchmark-paper

# Create mamba environment
mamba env create -f environment.yaml

# Activate environment
mamba activate q100-smk

# Verify installation
snakemake --version  # Should be ≥8.0
snakefmt --version
```

### Daily Workflow

```bash
# Activate environment before starting work
mamba activate q100-smk

# Update environment if dependencies changed
mamba env update -f environment.yaml
```

## Pipeline Development

### Running the Pipeline

```bash
# Dry-run to validate workflow
snakemake -n

# Dry-run with quiet output (less verbose)
snakemake -n --quiet

# Run pipeline locally with conda
snakemake --cores 4 --sdm conda

# Run with containers for reproducibility
snakemake --cores 4 --sdm conda apptainer

# Override default resources
snakemake --cores 4 --sdm conda \
    --set-resources generate_sv_len:mem_mb=8192

# Force re-run specific rule
snakemake --cores 4 --sdm conda \
    --forcerun generate_sv_len
```

### Testing and Validation

```bash
# Lint workflow (check best practices)
snakemake --lint

# Format Snakemake files
snakefmt workflow/

# Check formatting without modifying
snakefmt --check workflow/

# Validate configuration schema
snakemake --configfile config/config.yaml -n

# Run smoke test (via Makefile)
make dry-run
make lint
```

### Common Issues

**Conda solver taking too long:**

```bash
# Use mamba instead of conda
snakemake --cores 4 --sdm mamba
```

**Logs not visible:**

```bash
# Check log directory exists
ls -la logs/

# View specific rule log
cat logs/vcf_processing/sv_len.log
```

**Output files not generated:**

```bash
# Check if input files exist
snakemake --summary

# Dry-run with full DAG
snakemake -n --printshellcmds
```

## Git Workflow

### Branch Management

```bash
# Always work on feature branches
git checkout main
git pull origin main

# Create feature branch (use conventional naming)
git checkout -b feat/add-rtg-stats

# Make changes, commit frequently
git add workflow/rules/vcf_processing.smk
git commit -m "feat: add RTG vcfstats rule for historical benchmarks"

# Push branch
git push -u origin feat/add-rtg-stats

# Open pull request on GitHub
```

### Commit Conventions

Use [Conventional Commits](https://www.conventionalcommits.org/) format:

**Types:**

- `feat:` - New feature or analysis
- `fix:` - Bug fix
- `docs:` - Documentation updates
- `chore:` - Maintenance (dependencies, configs)
- `refactor:` - Code restructuring
- `test:` - Adding/updating tests

**Examples:**

```bash
# New feature
git commit -m "feat: add bedtools stratification intersect rule"

# Bug fix
git commit -m "fix: correct SVLEN extraction in bcftools query"

# Documentation
git commit -m "docs: add bedtools usage examples to README"

# Chore
git commit -m "chore: pin bcftools to version 1.19"

# Multi-line commit with body
git commit -m "feat: add historical benchmark downloads

- Added v4.2.1 GRCh38 VCF/BED downloads
- Added CMRG v1.00 SV benchmark downloads
- Configured SHA256 checksums for validation
"
```

### Before Committing

```bash
# Run linting and formatting
snakemake --lint
snakefmt workflow/

# Dry-run to validate workflow
snakemake -n --quiet

# Check what will be committed
git status
git diff

# Stage and commit
git add <files>
git commit -m "type: description"
```

## Tool Documentation

### bcftools

**Documentation:** [https://samtools.github.io/bcftools/](https://samtools.github.io/bcftools/)

**Common commands:**

```bash
# Query variant fields
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' input.vcf.gz

# Filter by region
bcftools view -r chr1:1000-2000 input.vcf.gz

# Statistics
bcftools stats input.vcf.gz

# Index VCF
bcftools index input.vcf.gz
```

### bedtools

**Documentation:** [https://bedtools.readthedocs.io/](https://bedtools.readthedocs.io/)

**Common commands:**

```bash
# Intersect regions
bedtools intersect -a regions.bed -b features.bed

# Calculate coverage
bedtools coverage -a regions.bed -b features.bed

# Genome-wide coverage
bedtools genomecov -i regions.bed -g genome.fa.fai

# Merge overlapping regions
bedtools merge -i regions.bed
```

### seqtk

**Documentation:** [https://github.com/lh3/seqtk](https://github.com/lh3/seqtk)

**Common commands:**

```bash
# Compute sequence composition
seqtk comp genome.fa

# Extract subsequence
seqtk subseq genome.fa regions.bed

# Sample reads
seqtk sample -s100 reads.fq 10000 > subset.fq
```

### Snakemake

**Documentation:** [https://snakemake.readthedocs.io/](https://snakemake.readthedocs.io/)

**Wrapper Catalog:** [https://snakemake-wrappers.readthedocs.io/](https://snakemake-wrappers.readthedocs.io/)

**Key concepts:**

- **Rules:** Define workflow steps
- **Wildcards:** Generalize patterns (`{sample}`, `{ref}`)
- **Config:** External configuration via YAML
- **Conda/Container:** Isolated environments per rule
- **DAG:** Directed acyclic graph of dependencies

## Code Review Checklist

Before requesting review:

- [ ] All rules have `log:`, `threads:`, `resources:`, `conda:` directives
- [ ] Outputs use `ensure()` for validation
- [ ] Download rules have `retries: 3` and `ancient()`
- [ ] Conda envs exclude `defaults` channel
- [ ] `snakemake --lint` passes
- [ ] `snakefmt --check workflow/` passes
- [ ] `snakemake -n` completes without errors
- [ ] Commit messages follow Conventional Commits
- [ ] README updated if adding new features

## Additional Resources

- **Snakemake Tutorial:** https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html
- **Best Practices:** https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html
- **Bioconda:** https://bioconda.github.io/
- **Conda-forge:** https://conda-forge.org/
