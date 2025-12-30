# Plan: Switch to bcftools query for Variant Counts

This plan outlines the steps to modify the Snakemake pipeline to generate variant count tables using `bcftools query` instead of relying solely on `rtg vcfstats`. This approach allows for more flexible downstream analysis in R/Quarto.

## 1. Config Updates

### `config/config.yaml`
-   **Action**: Add a new top-level section `stratifications`.
-   **Purpose**: Define stratification regions (difficult contexts) dynamically instead of hardcoding them in `common.smk`.
-   **Structure**:
    ```yaml
    stratifications:
      TR:
        path_template: "LowComplexity/{ref}_AllTandemRepeats.bed.gz"
        description: "Tandem Repeats"
      HP:
        path_template: "LowComplexity/{ref}_Homopolymers.bed.gz"
        description: "Homopolymers"
      SD:
        path_template: "SegDup/{ref}_SegDups.bed.gz"
        description: "Segmental Duplications"
      # Add other stratifications as needed (e.g., GC content, Mappability)
    ```
-   **Note**: The `path_template` uses `{ref}` which will be substituted (e.g., `GRCh38` or `hg19`) based on the benchmark's reference.

## 2. Snakemake Rule Updates

### `workflow/rules/common.smk`
-   **Action**: Update logic to read stratifications from `config['stratifications']`.
-   **Details**:
    -   Remove hardcoded `CONTEXT_PATHS`.
    -   Add helper functions:
        -   `get_stratification_path(strat_name, ref_version)`: Resolves the URL/path for a given stratification and reference.
        -   `get_stratification_names()`: Returns list of keys from `config['stratifications']`.

### `workflow/rules/vcf_processing.smk`

#### Update `download_stratification`
-   **Action**: Modify to use the new config structure.
-   **Input**: Uses the URL templates from `config['stratifications']`.
-   **Output**: `resources/stratifications/{ref}/{strat_name}.bed.gz`.

#### New Rule: `annotate_benchmark_vcf`
-   **Purpose**: Create a single VCF for each benchmark that contains all necessary context information as INFO tags.
-   **Input**:
    -   Benchmark VCF: `resources/benchmarksets/{benchmark}/{benchmark}.vcf.gz`
    -   Benchmark BED: `resources/benchmarksets/{benchmark}/{benchmark}.bed.gz`
    -   Stratification BEDs: `expand("resources/stratifications/{ref}/{strat}.bed.gz", ...)`
    -   Exclusion BEDs (for v5q benchmarks): `results/exclusions/{benchmark}_exclusions.bed.gz` (if applicable).
-   **Output**: `results/annotated_vcfs/{benchmark}/{benchmark}.annotated.vcf.gz` (and index).
-   **Tools**: `bcftools annotate`.
-   **Steps**:
    1.  **Header Preparation**: Create a header file defining new INFO tags:
        -   `INFO/BMKREGIONS`: Flag, "Variant falls within benchmark regions".
        -   `INFO/STRAT_{NAME}`: Flag, "Variant overlaps {Description} (>20%)".
        -   `INFO/EXCL_{NAME}`: Flag, "Variant overlaps exclusion {Name}".
    2.  **Annotation**:
        -   Annotate with Benchmark BED -> `BMKREGIONS`.
        -   Annotate with Stratification BEDs -> `STRAT_{NAME}`. **Crucial**: Use `--min-overlap 0.20` for stratifications to match GA4GH standards for difficult regions.
        -   Annotate with Exclusion BEDs -> `EXCL_{NAME}` (exact overlap usually sufficient, or define policy).
    3.  **Output**: Write compressed VCF.

#### New Rule: `generate_variant_table`
-   **Purpose**: Extract variant-level data into a tabular format for analysis.
-   **Input**: `results/annotated_vcfs/{benchmark}/{benchmark}.annotated.vcf.gz`.
-   **Output**: `results/var_tables/{benchmark}/{benchmark}.tsv.gz`.
-   **Tools**: `bcftools norm`, `bcftools query`, python script.
-   **Steps**:
    1.  **Split Multiallelics**: Pipe input to `bcftools norm -m -any`. This ensures one line per allele, simplifying length calculations.
    2.  **Query**: Run `bcftools query` to extract:
        -   CHROM, POS, REF, ALT, QUAL, FILTER
        -   GT (Genotype)
        -   INFO tags: `BMKREGIONS`, `STRAT_*`, `EXCL_*`
    3.  **Format**: Pipe output to `workflow/scripts/format_variant_table.py`.
    4.  **Compress**: Output as `.tsv.gz`.

## 3. Helper Scripts

### `workflow/scripts/format_variant_table.py`
-   **Purpose**: Process the raw stream from `bcftools query` and calculate derived fields.
-   **Logic**:
    -   Read stdin line by line.
    -   **Calculate Lengths**:
        -   `ref_len = len(REF)`
        -   `alt_len = len(ALT)`
        -   `ilen = alt_len - ref_len` (Insertion > 0, Deletion < 0, SNP = 0).
        -   `abs_ilen = abs(ilen)`
    -   **Determine Type**:
        -   SNP: `ilen == 0` and `ref_len == 1` and `alt_len == 1`.
        -   INDEL: `ilen != 0`.
        -   COMPLEX: Other cases (MNP, etc., if not normalized).
    -   **Format Boolean Flags**: Convert `.` (missing) to `0` or `False`, and present tags to `1` or `True` for the INFO columns.
    -   **Write Output**: Tab-separated values with a header.

## 4. Documentation

### `README.qmd`
-   Update the "Outputs" section to describe:
    -   `results/annotated_vcfs/`: Intermediate VCFs with context annotations.
    -   `results/var_tables/`: Final TSVs used for manuscript figures.
-   Note that `rtg_vcfstats` (if kept) provides summary stats, but `var_tables` are the source of truth for detailed breakdowns.

## 5. Implementation Details

### TSV Columns
The final TSV should contain at least:
-   `chrom`: Chromosome
-   `pos`: Position (1-based)
-   `ref`: Reference allele
-   `alt`: Alternate allele
-   `gt`: Genotype (e.g., 0/1, 1/1)
-   `type`: SNP, INDEL, COMPLEX
-   `ilen`: Insertion/Deletion length (negative for deletions)
-   `abs_ilen`: Absolute indel length
-   `in_benchmark`: Boolean (from `BMKREGIONS`)
-   `strat_tr`: Boolean (Tandem Repeats)
-   `strat_hp`: Boolean (Homopolymers)
-   `strat_sd`: Boolean (Seg Dups)
-   ... (other stratifications)
-   `excl_all`: Boolean (Any exclusion)
-   ... (specific exclusions if needed)

### Naming Conventions
-   **INFO Tags**:
    -   Stratifications: `STRAT_{UPPERCASE_NAME}` (e.g., `STRAT_TR`).
    -   Exclusions: `EXCL_{UPPERCASE_NAME}`.
    -   Benchmark Regions: `BMKREGIONS`.
