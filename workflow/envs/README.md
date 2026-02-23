# Conda Environment Structure

## Active Environments (5)

### biotools.yaml
**Purpose:** Core CLI bioinformatics tools
**Key packages:** bcftools=1.22, rtg-tools=3.13
**Used by:** VCF processing, benchmark comparisons
**Replaces:** bcftools.yaml + rtg-tools.yaml (consolidated 2026-02-23)

### python-biotools.yaml
**Purpose:** Python data processing + genomic interval tools
**Key packages:** python=3.11, pandas, pyarrow, bcftools=1.22, bedtools=2.31.1, tabix
**Used by:** Annotation, genomic context analysis, exclusions pipeline
**Replaces:** python.yaml + bedtools.yaml (consolidated 2026-02-23)

### samtools.yaml
**Purpose:** Sequence handling and manipulation
**Key packages:** samtools=1.22, seqkit=2.12.0
**Used by:** Reference processing, chr8 synteny pipeline
**Status:** Kept separate (lightweight, clear separation of concerns)

### downloads.yaml
**Purpose:** File retrieval from remote sources
**Key packages:** wget
**Used by:** Download rules for references, stratifications, benchmarks
**Status:** Kept separate (minimal environment, fast builds)

### plotsr.yaml
**Purpose:** Chr8 synteny visualization pipeline
**Key packages:** minimap2, samtools, syri, plotsr, matplotlib, **pandas<2.0**
**Used by:** Chr8 synteny analysis and figure generation
**Status:** **MANDATORY ISOLATION** - Must remain separate due to pandas<2.0 constraint (SyRI Cython bug with pandas 2.0+)

### truvari.yaml
**Purpose:** Variant analysis and comparison
**Key packages:** Truvari==5.4.0 (pip), bcftools=1.20, bedtools, pandas, pyarrow=14.0
**Used by:** Variant table generation, variant counting, benchmark comparisons
**Status:** Kept separate due to bcftools version conflict (uses 1.20 vs 1.22 in biotools)

---

## Deprecated Environments

See `deprecated/` for historical environments removed during consolidation:
- **bcftools.yaml** - Merged into biotools.yaml
- **python.yaml** - Merged into python-biotools.yaml
- **bedtools.yaml** - Merged into python-biotools.yaml
- **rtg-tools.yaml** - Merged into biotools.yaml

**Date deprecated:** 2026-02-23
**Consolidation impact:** 8 → 6 environments, ~800 MB disk space savings

---

## Consolidation Benefits

1. **Reduced duplication:** Eliminated redundant package installations (bcftools in 4 envs → 3 envs, Python/pandas in 4 envs → 2 envs)
2. **Faster builds:** Fewer environments to solve and build in CI/clean installs
3. **Easier maintenance:** Updating shared dependencies (bcftools, pandas) requires fewer environment files
4. **Preserved isolation:** Critical constraints (plotsr pandas<2.0, truvari bcftools=1.20) remain isolated

---

## Testing Protocol

Before modifying environments, run full test suite:

1. **Environment build test:**
   ```bash
   mamba env create -f workflow/envs/<env-name>.yaml -n test-<env> --dry-run
   ```

2. **Dry-run test:**
   ```bash
   make dry-run
   ```

3. **Integration test** (single benchmark):
   ```bash
   snakemake --cores 4 --use-conda \
     results/variant_tables/v5.0q_GRCh38_smvar/variants.parquet
   ```

4. **Full pipeline validation:**
   ```bash
   make clean && make run
   ```

---

## Future Consolidation Opportunities

**Potential:** Investigate merging truvari.yaml into python-biotools.yaml
- **Blocker:** bcftools version conflict (1.20 vs 1.22)
- **Risk:** Truvari may have undocumented bcftools=1.20 dependency
- **Testing required:**
  1. Test Truvari with bcftools=1.22
  2. Verify variant classification unchanged
  3. Compare outputs byte-for-byte with baseline

**Not recommended:**
- Merging plotsr.yaml (pandas<2.0 is hard constraint)
- Merging downloads.yaml (no benefit, adds complexity)
