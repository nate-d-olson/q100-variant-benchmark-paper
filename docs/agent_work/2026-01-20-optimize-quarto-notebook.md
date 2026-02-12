# Optimize Quarto Notebook by Moving Data Processing to Snakemake

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Speed up benchmarkset-characterization.qmd by moving heavy file parsing and data transformation from R to Snakemake, reducing R memory usage and computation time.

**Architecture:** Create a new Snakemake rule that generates pre-processed variant tables with cleaned data (variant type corrections, variant size annotations, chromosome formatting) that R can directly read without additional processing.

**Tech Stack:** Python (pandas), Snakemake, bcftools/existing pipeline infrastructure

---

## Task 1: Create Python Script for Variant Table Processing

**Files:**

- Create: `workflow/scripts/process_variant_table.py`

**Step 1: Write failing test for variant table processing**

```python
# tests/test_process_variant_table.py
import pandas as pd
import pytest
from pathlib import Path

def test_tidy_smvar_filters_benchmark_regions():
    """Test that tidy_smvar removes variants outside benchmark regions."""
    df = pd.DataFrame({
        'chrom': ['chr1', 'chr1', 'chr2'],
        'pos': [100, 200, 300],
        'end': [101, 201, 301],
        'gt': ['0/1', '1/1', '0/1'],
        'vkx': ['var1', 'var2', 'var3'],
        'var_type': ['SNP', 'INDEL', 'SNP'],
        'len_ref': [1, 5, 1],
        'len_alt': [1, 1, 1],
        'strat_ids': ['.', '.', '.'],
        'region_ids': ['BMKREGIONS', 'OTHER', 'BMKREGIONS']
    })

    from workflow.scripts.process_variant_table import tidy_smvar

    result = tidy_smvar(df)

    assert len(result) == 2
    assert all(result['region_ids'].str.startswith('BMKREGIONS'))

def test_tidy_smvar_filters_size_lt50():
    """Test that tidy_smvar keeps only variants < 50bp."""
    df = pd.DataFrame({
        'chrom': ['chr1', 'chr1', 'chr2'],
        'pos': [100, 200, 300],
        'end': [101, 201, 301],
        'gt': ['0/1', '1/1', '0/1'],
        'vkx': ['var1', 'var2', 'var3'],
        'var_type': ['SNP', 'INDEL', 'INDEL'],
        'len_ref': [1, 5, 60],
        'len_alt': [1, 1, 1],
        'strat_ids': ['.', '.', '.'],
        'region_ids': ['BMKREGIONS', 'BMKREGIONS', 'BMKREGIONS']
    })

    from workflow.scripts.process_variant_table import tidy_smvar

    result = tidy_smvar(df)

    assert len(result) == 2
    assert all((result['len_ref'] < 50) & (result['len_alt'] < 50))

def test_tidy_smvar_corrects_variant_types():
    """Test that tidy_smvar corrects OTHER and OVERLAP variant types."""
    df = pd.DataFrame({
        'chrom': ['chr1', 'chr1', 'chr2', 'chr3'],
        'pos': [100, 200, 300, 400],
        'end': [101, 201, 301, 401],
        'gt': ['0/1', '1/1', '0/1', '1/1'],
        'vkx': ['var1', 'var2', 'var3', 'var4'],
        'var_type': ['OTHER', 'OVERLAP', 'OVERLAP', 'OVERLAP'],
        'len_ref': [5, 1, 5, 5],
        'len_alt': [1, 1, 1, 5],
        'strat_ids': ['.', '.', '.', '.'],
        'region_ids': ['BMKREGIONS', 'BMKREGIONS', 'BMKREGIONS', 'BMKREGIONS']
    })

    from workflow.scripts.process_variant_table import tidy_smvar

    result = tidy_smvar(df)

    # OTHER with len_ref>1, len_alt=1 should be INDEL
    assert result.loc[result['vkx'] == 'var1', 'var_type'].iloc[0] == 'INDEL'
    # OVERLAP with len_ref=1, len_alt=1 should be SNP
    assert result.loc[result['vkx'] == 'var2', 'var_type'].iloc[0] == 'SNP'
    # OVERLAP with len_ref>1, len_alt=1 should be INDEL
    assert result.loc[result['vkx'] == 'var3', 'var_type'].iloc[0] == 'INDEL'
    # OVERLAP with len_ref>1, len_alt>1 should be COMPLEX
    assert result.loc[result['vkx'] == 'var4', 'var_type'].iloc[0] == 'COMPLEX'

def test_tidy_smvar_adds_var_size():
    """Test that tidy_smvar adds var_size column."""
    df = pd.DataFrame({
        'chrom': ['chr1', 'chr1'],
        'pos': [100, 200],
        'end': [101, 201],
        'gt': ['0/1', '1/1'],
        'vkx': ['var1', 'var2'],
        'var_type': ['SNP', 'INDEL'],
        'len_ref': [1, 5],
        'len_alt': [1, 2],
        'strat_ids': ['.', '.'],
        'region_ids': ['BMKREGIONS', 'BMKREGIONS']
    })

    from workflow.scripts.process_variant_table import tidy_smvar

    result = tidy_smvar(df)

    assert 'var_size' in result.columns
    assert result.loc[result['vkx'] == 'var1', 'var_size'].iloc[0] == 0
    assert result.loc[result['vkx'] == 'var2', 'var_size'].iloc[0] == -3

def test_tidy_stvar_filters_benchmark_regions():
    """Test that tidy_stvar removes variants outside benchmark regions."""
    df = pd.DataFrame({
        'chrom': ['chr1', 'chr1'],
        'pos': [100, 200],
        'end': [200, 300],
        'gt': ['0/1', '1/1'],
        'vkx': ['var1', 'var2'],
        'var_type': ['DEL', 'INS'],
        'len_ref': [100, 1],
        'len_alt': [1, 100],
        'strat_ids': ['.', '.'],
        'region_ids': ['BMKREGIONS', 'OTHER'],
        'SVTYPE': ['DEL', 'INS'],
        'SVLEN': [99, 99]
    })

    from workflow.scripts.process_variant_table import tidy_stvar

    result = tidy_stvar(df)

    assert len(result) == 1
    assert all(result['region_ids'].str.startswith('BMKREGIONS'))

def test_tidy_stvar_filters_size_gt50():
    """Test that tidy_stvar keeps only variants >= 50bp."""
    df = pd.DataFrame({
        'chrom': ['chr1', 'chr1', 'chr2'],
        'pos': [100, 200, 300],
        'end': [200, 250, 400],
        'gt': ['0/1', '1/1', '0/1'],
        'vkx': ['var1', 'var2', 'var3'],
        'var_type': ['DEL', 'INS', 'DEL'],
        'len_ref': [100, 1, 10],
        'len_alt': [1, 100, 1],
        'strat_ids': ['.', '.', '.'],
        'region_ids': ['BMKREGIONS', 'BMKREGIONS', 'BMKREGIONS'],
        'SVTYPE': ['DEL', 'INS', 'DEL'],
        'SVLEN': [99, 99, 9]
    })

    from workflow.scripts.process_variant_table import tidy_stvar

    result = tidy_stvar(df)

    assert len(result) == 2
    assert all((result['len_ref'] > 49) | (result['len_alt'] > 49))

def test_tidy_stvar_sets_var_type_from_svtype():
    """Test that tidy_stvar uses SVTYPE for var_type."""
    df = pd.DataFrame({
        'chrom': ['chr1', 'chr1'],
        'pos': [100, 200],
        'end': [200, 300],
        'gt': ['0/1', '1/1'],
        'vkx': ['var1', 'var2'],
        'var_type': ['OTHER', 'COMPLEX'],
        'len_ref': [100, 1],
        'len_alt': [1, 100],
        'strat_ids': ['.', '.'],
        'region_ids': ['BMKREGIONS', 'BMKREGIONS'],
        'SVTYPE': ['DEL', 'INS'],
        'SVLEN': [99, 99]
    })

    from workflow.scripts.process_variant_table import tidy_stvar

    result = tidy_stvar(df)

    assert result.loc[result['vkx'] == 'var1', 'var_type'].iloc[0] == 'DEL'
    assert result.loc[result['vkx'] == 'var2', 'var_type'].iloc[0] == 'INS'

def test_tidy_stvar_calculates_var_size():
    """Test that tidy_stvar calculates var_size from SVTYPE and SVLEN."""
    df = pd.DataFrame({
        'chrom': ['chr1', 'chr1'],
        'pos': [100, 200],
        'end': [200, 300],
        'gt': ['0/1', '1/1'],
        'vkx': ['var1', 'var2'],
        'var_type': ['DEL', 'INS'],
        'len_ref': [100, 1],
        'len_alt': [1, 100],
        'strat_ids': ['.', '.'],
        'region_ids': ['BMKREGIONS', 'BMKREGIONS'],
        'SVTYPE': ['DEL', 'INS'],
        'SVLEN': [99, 150]
    })

    from workflow.scripts.process_variant_table import tidy_stvar

    result = tidy_stvar(df)

    assert 'var_size' in result.columns
    assert result.loc[result['vkx'] == 'var1', 'var_size'].iloc[0] == -99
    assert result.loc[result['vkx'] == 'var2', 'var_size'].iloc[0] == 150

def test_format_chrom_adds_chr_prefix():
    """Test that format_chrom adds chr prefix for GRCh37."""
    df = pd.DataFrame({
        'chrom': ['1', '2', 'X', 'Y']
    })

    from workflow.scripts.process_variant_table import format_chrom

    result = format_chrom(df, is_grch37=True)

    assert all(result['chrom'].str.startswith('chr'))
    assert result['chrom'].tolist() == ['chr1', 'chr2', 'chrX', 'chrY']

def test_format_chrom_preserves_chr_prefix():
    """Test that format_chrom preserves chr prefix for GRCh38/CHM13."""
    df = pd.DataFrame({
        'chrom': ['chr1', 'chr2', 'chrX', 'chrY']
    })

    from workflow.scripts.process_variant_table import format_chrom

    result = format_chrom(df, is_grch37=False)

    assert result['chrom'].tolist() == ['chr1', 'chr2', 'chrX', 'chrY']

def test_format_chrom_creates_categorical():
    """Test that format_chrom creates categorical dtype with correct levels."""
    df = pd.DataFrame({
        'chrom': ['chr1', 'chr22', 'chrX', 'chrY']
    })

    from workflow.scripts.process_variant_table import format_chrom

    result = format_chrom(df, is_grch37=False)

    assert result['chrom'].dtype.name == 'category'
    expected_levels = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
    assert result['chrom'].cat.categories.tolist() == expected_levels

def test_select_output_columns():
    """Test that select_output_columns keeps only required columns."""
    df = pd.DataFrame({
        'chrom': ['chr1'],
        'pos': [100],
        'end': [101],
        'gt': ['0/1'],
        'vkx': ['var1'],
        'var_type': ['SNP'],
        'var_size': [0],
        'region_ids': ['BMKREGIONS'],
        'strat_ids': ['.'],
        'extra1': ['foo'],
        'extra2': ['bar']
    })

    from workflow.scripts.process_variant_table import select_output_columns

    result = select_output_columns(df)

    expected_cols = ['chrom', 'pos', 'end', 'gt', 'vkx', 'var_type',
                     'var_size', 'region_ids', 'strat_ids']
    assert result.columns.tolist() == expected_cols
```

**Step 2: Run tests to verify they fail**

Run: `pytest tests/test_process_variant_table.py -v`
Expected: FAIL with "ModuleNotFoundError: No module named 'workflow.scripts.process_variant_table'"

**Step 3: Write minimal implementation**

```python
# workflow/scripts/process_variant_table.py
"""
Process variant tables from Snakemake pipeline output.

This script replicates the tidy_smvar and tidy_stvar functions from
benchmarkset-characterization.qmd to pre-process variant tables in
the Snakemake pipeline.
"""

from pathlib import Path
from typing import Literal
import pandas as pd
import argparse
import sys


def tidy_smvar(df: pd.DataFrame) -> pd.DataFrame:
    """
    Clean and annotate small variant table.

    Filters:
    - Removes variants outside benchmark regions
    - Keeps only variants < 50bp

    Transformations:
    - Corrects var_type for OTHER and OVERLAP variants
    - Adds var_size column (len_alt - len_ref)

    Args:
        df: Variant table with standard columns

    Returns:
        Cleaned variant table
    """
    # Remove variants outside benchmark regions
    df = df[df['region_ids'].str.startswith('BMKREGIONS')].copy()

    # Keep only variants < 50bp
    df = df[(df['len_ref'] < 50) & (df['len_alt'] < 50)].copy()

    # Correct variant types
    def correct_var_type(row):
        var_type = row['var_type']
        len_ref = row['len_ref']
        len_alt = row['len_alt']

        # OTHER variants
        if var_type == 'OTHER':
            if (len_ref > 1 and len_alt == 1) or (len_ref == 1 and len_alt > 1):
                return 'INDEL'

        # OVERLAP variants
        elif var_type == 'OVERLAP':
            if len_ref == 1 and len_alt == 1:
                return 'SNP'
            elif (len_ref > 1 and len_alt == 1) or (len_ref == 1 and len_alt > 1):
                return 'INDEL'
            elif len_ref > 1 and len_alt > 1:
                return 'COMPLEX'

        return var_type

    df['var_type'] = df.apply(correct_var_type, axis=1)

    # Add variant size
    df['var_size'] = df['len_alt'] - df['len_ref']

    return df


def tidy_stvar(df: pd.DataFrame) -> pd.DataFrame:
    """
    Clean and annotate structural variant table.

    Filters:
    - Removes variants outside benchmark regions
    - Keeps only variants >= 50bp

    Transformations:
    - Sets var_type from SVTYPE
    - Calculates var_size from SVTYPE and SVLEN

    Args:
        df: Variant table with standard columns plus SVTYPE, SVLEN

    Returns:
        Cleaned variant table
    """
    # Remove variants outside benchmark regions
    df = df[df['region_ids'].str.startswith('BMKREGIONS')].copy()

    # Keep only variants >= 50bp
    df = df[(df['len_ref'] > 49) | (df['len_alt'] > 49)].copy()

    # Set var_type from SVTYPE
    df['var_type'] = df['SVTYPE']

    # Calculate var_size
    df['var_size'] = 0
    df.loc[df['SVTYPE'] == 'INS', 'var_size'] = df.loc[df['SVTYPE'] == 'INS', 'SVLEN']
    df.loc[df['SVTYPE'] == 'DEL', 'var_size'] = -df.loc[df['SVTYPE'] == 'DEL', 'SVLEN']

    return df


def format_chrom(df: pd.DataFrame, is_grch37: bool) -> pd.DataFrame:
    """
    Format chromosome names with chr prefix and convert to categorical.

    Args:
        df: DataFrame with 'chrom' column
        is_grch37: Whether this is GRCh37 (needs chr prefix added)

    Returns:
        DataFrame with formatted categorical chrom column
    """
    df = df.copy()

    # Add chr prefix for GRCh37
    if is_grch37:
        df['chrom'] = 'chr' + df['chrom'].astype(str)

    # Convert to categorical with proper ordering
    chrom_levels = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
    df['chrom'] = pd.Categorical(df['chrom'], categories=chrom_levels, ordered=True)

    return df


def select_output_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Select only the columns needed for downstream analysis.

    Args:
        df: Full variant table

    Returns:
        DataFrame with only required columns
    """
    output_cols = [
        'chrom', 'pos', 'end', 'gt', 'vkx',
        'var_type', 'var_size', 'region_ids', 'strat_ids'
    ]

    return df[output_cols]


def process_variant_table(
    input_path: Path,
    output_path: Path,
    variant_type: Literal['smvar', 'stvar'],
    is_grch37: bool
) -> None:
    """
    Process a variant table and write cleaned output.

    Args:
        input_path: Path to input TSV file
        output_path: Path to output TSV file
        variant_type: 'smvar' or 'stvar'
        is_grch37: Whether this is GRCh37 (for chrom formatting)
    """
    print(f"Reading {input_path}", file=sys.stderr)

    # Read TSV with pandas
    df = pd.read_csv(input_path, sep='\t', comment='#')

    # Clean column names (remove [#] prefix from bcftools query output)
    df.columns = df.columns.str.replace(r'^\[.*?\]', '', regex=True)

    print(f"Loaded {len(df)} variants", file=sys.stderr)

    # Apply appropriate cleaning function
    if variant_type == 'smvar':
        print("Tidying small variant table", file=sys.stderr)
        df = tidy_smvar(df)
    elif variant_type == 'stvar':
        print("Tidying structural variant table", file=sys.stderr)
        df = tidy_stvar(df)
    else:
        raise ValueError(f"Unknown variant type: {variant_type}")

    print(f"After filtering: {len(df)} variants", file=sys.stderr)

    # Format chromosomes
    df = format_chrom(df, is_grch37=is_grch37)

    # Select output columns
    df = select_output_columns(df)

    # Write output
    print(f"Writing {output_path}", file=sys.stderr)
    df.to_csv(output_path, sep='\t', index=False)

    print(f"Done. Wrote {len(df)} variants", file=sys.stderr)


def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description='Process variant tables from Snakemake pipeline'
    )
    parser.add_argument(
        '--input', '-i',
        type=Path,
        required=True,
        help='Input TSV file'
    )
    parser.add_argument(
        '--output', '-o',
        type=Path,
        required=True,
        help='Output TSV file'
    )
    parser.add_argument(
        '--type', '-t',
        choices=['smvar', 'stvar'],
        required=True,
        help='Variant type'
    )
    parser.add_argument(
        '--grch37',
        action='store_true',
        help='Input is GRCh37 (needs chr prefix added)'
    )

    args = parser.parse_args()

    process_variant_table(
        input_path=args.input,
        output_path=args.output,
        variant_type=args.type,
        is_grch37=args.grch37
    )


if __name__ == '__main__':
    main()
```

**Step 4: Run tests to verify they pass**

Run: `pytest tests/test_process_variant_table.py -v`
Expected: All tests PASS

**Step 5: Commit**

```bash
git add workflow/scripts/process_variant_table.py tests/test_process_variant_table.py
git commit -m "feat: add variant table processing script

- Implements tidy_smvar and tidy_stvar from R
- Filters variants by benchmark regions and size
- Corrects variant types for OTHER and OVERLAP
- Adds var_size annotation
- Formats chromosome names as categorical
- Includes comprehensive test suite"
```

---

## Task 2: Add Snakemake Rule for Processing Variant Tables

**Files:**

- Modify: `workflow/rules/var_tables.smk`
- Create: `workflow/envs/pandas.yaml`

**Step 1: Write test for processed variant table output**

Create manual test to verify rule produces expected output:

```bash
# tests/test_process_variant_table_rule.sh
#!/bin/bash
# Test that the process_variant_table rule produces valid output

set -e

# Run on a single benchmark
snakemake -n results/variant_tables_processed/v5q_GRCh38_smvar/variants_processed.tsv

# Check that rule is defined
if [ $? -eq 0 ]; then
    echo "PASS: Rule is defined"
else
    echo "FAIL: Rule not found"
    exit 1
fi
```

**Step 2: Run test to verify it fails**

Run: `bash tests/test_process_variant_table_rule.sh`
Expected: FAIL with rule not found

**Step 3: Implement the Snakemake rule**

First create pandas environment:

```yaml
# workflow/envs/pandas.yaml
channels:
  - conda-forge
  - bioconda
dependencies:
  - python>=3.10
  - pandas>=2.0
  - pytest>=7.0
```

Then add rule to var_tables.smk:

```python
# Add to workflow/rules/var_tables.smk

rule process_variant_table:
    """
    Process variant table: filter, correct types, add annotations.

    This replaces the R processing in benchmarkset-characterization.qmd.
    """
    input:
        tsv="results/variant_tables/{benchmark}/variants.tsv",
    output:
        tsv=ensure(
            "results/variant_tables_processed/{benchmark}/variants_processed.tsv",
            non_empty=True
        ),
    params:
        var_type=lambda w: "stvar" if w.benchmark.startswith("v06_") or "stvar" in w.benchmark else "smvar",
        is_grch37=lambda w: "--grch37" if "GRCh37" in w.benchmark else "",
    log:
        "logs/process_variant_table/{benchmark}.log",
    conda:
        "../envs/pandas.yaml"
    shell:
        """
        python workflow/scripts/process_variant_table.py \
            --input {input.tsv} \
            --output {output.tsv} \
            --type {params.var_type} \
            {params.is_grch37} 2> {log}
        """
```

**Step 4: Update common.smk to include processed tables in outputs**

```python
# Modify workflow/rules/common.smk

def get_var_table_inputs(wildcards) -> List[str]:
    """
    Generate list of variant table files for all benchmarks.

    Returns:
        List of file paths for variant table outputs
    """
    outputs = []
    for benchmark in config["benchmarksets"]:
        # Raw tables
        outputs.append(f"results/variant_tables/{benchmark}/variants.tsv")
        # Processed tables
        outputs.append(f"results/variant_tables_processed/{benchmark}/variants_processed.tsv")
    return outputs
```

**Step 5: Test the rule**

Run: `snakemake -n results/variant_tables_processed/v5q_GRCh38_smvar/variants_processed.tsv`
Expected: Rule executes successfully and produces valid TSV

**Step 6: Commit**

```bash
git add workflow/rules/var_tables.smk workflow/envs/pandas.yaml workflow/rules/common.smk
git commit -m "feat: add Snakemake rule to process variant tables

- Adds process_variant_table rule
- Creates pandas conda environment
- Integrates processed tables into pipeline outputs
- Replaces R processing with Python preprocessing"
```

---

## Task 3: Update Quarto Notebook to Use Processed Tables

**Files:**

- Modify: `analysis/benchmarkset-characterization.qmd`

**Step 1: Write test to verify notebook runs faster**

```bash
# tests/test_notebook_performance.sh
#!/bin/bash
# Test that updated notebook runs successfully

set -e

# Time the notebook rendering
echo "Testing notebook with processed tables..."
time quarto render analysis/benchmarkset-characterization.qmd

if [ $? -eq 0 ]; then
    echo "PASS: Notebook rendered successfully"
else
    echo "FAIL: Notebook rendering failed"
    exit 1
fi
```

**Step 2: Run test to establish baseline**

Run: `bash tests/test_notebook_performance.sh`
Expected: Records baseline time (likely slow with current implementation)

**Step 3: Simplify the R code to read processed tables**

Replace the complex reading logic (lines 54-320) with simple reads:

```r
# analysis/benchmarkset-characterization.qmd

# Replace lines 54-320 with:

```{r}
#| warning: false
# Read pre-processed variant tables
library(tidyverse)
library(here)
library(vroom)

# Get processed table files
processed_var_tbl_files <- list.files(
    path = here("results/variant_tables_processed"),
    pattern = "variants_processed.tsv",
    recursive = TRUE,
    full.names = TRUE
)

# Derive labels from file paths
(names(processed_var_tbl_files) <- processed_var_tbl_files %>%
    str_extract("(?<=variant_tables_processed/).*(?=/variants_processed.tsv)") %>%
    str_replace_all("^v421", "v4.2.1") %>%
    str_replace_all("^v06", "v0.6") %>%
    str_replace_all("^v5q", "v5.0q") %>%
    str_replace_all("grch", "GRCh") %>%
    str_replace_all("chm13", "CHM13v2.0"))

# Read all processed tables with simple vroom (no custom column logic needed)
bench_vars_tbls <- processed_var_tbl_files %>%
    enframe(name = "benchmarkset", value = "path") %>%
    mutate(
        data = map(path, ~ vroom(.x, delim = "\t", show_col_types = FALSE))
    ) %>%
    select(benchmarkset, data) %>%
    separate(
        col = benchmarkset,
        into = c("bench_version", "ref", "bench_type"),
        sep = "_",
        fill = "right"
    ) %>%
    unnest(data)
```

# Remove these functions (no longer needed)

# - tidy_smvar() (lines 56-87)

# - tidy_stvar() (lines 89-111)

# - get_bench_var_cols() (lines 113-234)

# - read_bench_file() (lines 237-285)

```

**Step 4: Update TODO comment about optimization**

```r
# Update line 29 from:
# - Optimize variant tables: data.table, database, arrow, etc.?

# To:
# - DONE: Optimized by moving processing to Snakemake pipeline
```

**Step 5: Test notebook runs successfully**

Run: `quarto render analysis/benchmarkset-characterization.qmd`
Expected: Renders successfully with same output but faster

**Step 6: Commit**

```bash
git add analysis/benchmarkset-characterization.qmd
git commit -m "refactor: simplify Quarto notebook using processed tables

- Replace custom R parsing with simple vroom reads
- Remove tidy_smvar, tidy_stvar, get_bench_var_cols, read_bench_file
- Processing now handled by Snakemake pipeline
- Reduces R memory usage and computation time"
```

---

## Task 4: Run Full Pipeline and Verify Outputs

**Files:**

- None (testing/validation task)

**Step 1: Clean existing outputs**

```bash
# Remove old variant table outputs to force regeneration
rm -rf results/variant_tables_processed/
```

**Step 2: Run pipeline to generate processed tables**

Run: `snakemake --cores 12 --use-conda`
Expected: Pipeline completes successfully, generates all processed tables

**Step 3: Verify processed table contents**

```bash
# Check that processed tables have expected structure
head -n 20 results/variant_tables_processed/v5q_GRCh38_smvar/variants_processed.tsv

# Verify column names
head -n 1 results/variant_tables_processed/v5q_GRCh38_smvar/variants_processed.tsv

# Expected columns: chrom, pos, end, gt, vkx, var_type, var_size, region_ids, strat_ids
```

Expected: Tables have correct columns and data

**Step 4: Render Quarto notebook**

Run: `quarto render analysis/benchmarkset-characterization.qmd`
Expected: Renders successfully with identical output to original

**Step 5: Compare performance**

```bash
# Time the rendering
time quarto render analysis/benchmarkset-characterization.qmd
```

Expected: Significantly faster than baseline (check against baseline from Task 3)

**Step 6: Commit verification script**

```bash
# Create verification script for future use
cat > scripts/verify_processed_tables.sh << 'EOF'
#!/bin/bash
# Verify that processed variant tables are correctly formatted

set -e

echo "Verifying processed variant tables..."

for table in results/variant_tables_processed/*/variants_processed.tsv; do
    echo "Checking $table"

    # Check file exists and is non-empty
    if [ ! -s "$table" ]; then
        echo "ERROR: $table is empty or missing"
        exit 1
    fi

    # Check for expected columns
    header=$(head -n 1 "$table")
    expected="chrom pos end gt vkx var_type var_size region_ids strat_ids"

    if [ "$header" != "$expected" ]; then
        echo "ERROR: $table has unexpected columns"
        echo "Expected: $expected"
        echo "Got: $header"
        exit 1
    fi

    echo "  ✓ Valid"
done

echo "All processed tables verified successfully!"
EOF

chmod +x scripts/verify_processed_tables.sh
```

Run: `bash scripts/verify_processed_tables.sh`
Expected: All tables pass verification

```bash
git add scripts/verify_processed_tables.sh
git commit -m "test: add verification script for processed tables

- Checks that processed tables exist and are non-empty
- Verifies column names match expected format
- Can be run before rendering Quarto notebooks"
```

---

## Task 5: Update Documentation

**Files:**

- Modify: `README.md` or workflow documentation

**Step 1: Document the optimization**

Add section to README explaining the workflow:

```markdown
# README.md (add new section)

## Variant Table Processing

The pipeline generates two types of variant tables:

1. **Raw tables** (`results/variant_tables/{benchmark}/variants.tsv`)
   - Direct output from `bcftools query`
   - Contains all INFO fields from VCF

2. **Processed tables** (`results/variant_tables_processed/{benchmark}/variants_processed.tsv`)
   - Filtered to benchmark regions only
   - Size-filtered (< 50bp for smvar, >= 50bp for stvar)
   - Variant types corrected (OTHER → INDEL, OVERLAP → SNP/INDEL/COMPLEX)
   - var_size annotation added
   - Chromosome names formatted as categorical
   - Reduced to essential columns for analysis

The processed tables are used directly by the Quarto analysis notebooks,
eliminating the need for data processing in R and improving performance.

### Processing Logic

**Small variants (smvar):**
- Filters: in benchmark regions, both REF and ALT < 50bp
- Type corrections: OTHER and OVERLAP variants reclassified as SNP/INDEL/COMPLEX
- var_size: ALT length - REF length

**Structural variants (stvar):**
- Filters: in benchmark regions, REF or ALT >= 50bp
- var_type: Set from SVTYPE field
- var_size: SVLEN for INS (positive), -SVLEN for DEL (negative)
```

**Step 2: Update workflow diagram (if exists)**

If there's a workflow diagram, add the new processing step:

```
variant_tables → process_variant_table → variant_tables_processed → Quarto analysis
```

**Step 3: Commit documentation**

```bash
git add README.md
git commit -m "docs: document variant table processing optimization

- Explain raw vs processed tables
- Document filtering and transformation logic
- Update workflow description"
```

---

## Task 6: Optional - Add Performance Benchmarking

**Files:**

- Create: `tests/benchmark_notebook_performance.sh`

**Step 1: Create benchmarking script**

```bash
# tests/benchmark_notebook_performance.sh
#!/bin/bash
# Benchmark Quarto notebook rendering performance

set -e

ITERATIONS=3
OUTPUT_FILE="benchmark_results.txt"

echo "Benchmarking Quarto notebook rendering..." | tee "$OUTPUT_FILE"
echo "Iterations: $ITERATIONS" | tee -a "$OUTPUT_FILE"
echo "" | tee -a "$OUTPUT_FILE"

total_time=0

for i in $(seq 1 $ITERATIONS); do
    echo "Run $i/$ITERATIONS..." | tee -a "$OUTPUT_FILE"

    # Clear Quarto cache
    rm -rf analysis/benchmarkset-characterization_files/

    # Time the rendering
    start=$(date +%s)
    quarto render analysis/benchmarkset-characterization.qmd --quiet
    end=$(date +%s)

    elapsed=$((end - start))
    total_time=$((total_time + elapsed))

    echo "  Time: ${elapsed}s" | tee -a "$OUTPUT_FILE"
done

avg_time=$((total_time / ITERATIONS))

echo "" | tee -a "$OUTPUT_FILE"
echo "Average time: ${avg_time}s" | tee -a "$OUTPUT_FILE"
echo "Total time: ${total_time}s" | tee -a "$OUTPUT_FILE"
```

**Step 2: Run baseline benchmark (before optimization)**

```bash
# Checkout previous commit before optimization
git stash
git checkout HEAD~3  # Or appropriate commit before changes

chmod +x tests/benchmark_notebook_performance.sh
bash tests/benchmark_notebook_performance.sh

# Save baseline results
mv benchmark_results.txt benchmark_results_baseline.txt

# Return to current branch
git checkout -
git stash pop
```

**Step 3: Run optimized benchmark**

```bash
bash tests/benchmark_notebook_performance.sh

# Save optimized results
mv benchmark_results.txt benchmark_results_optimized.txt
```

**Step 4: Compare results**

```bash
echo "=== Performance Comparison ==="
echo "Baseline:"
cat benchmark_results_baseline.txt | grep "Average time"
echo ""
echo "Optimized:"
cat benchmark_results_optimized.txt | grep "Average time"
echo ""

# Calculate speedup
baseline=$(cat benchmark_results_baseline.txt | grep "Average time" | awk '{print $3}' | sed 's/s//')
optimized=$(cat benchmark_results_optimized.txt | grep "Average time" | awk '{print $3}' | sed 's/s//')

speedup=$(echo "scale=2; $baseline / $optimized" | bc)
echo "Speedup: ${speedup}x"
```

**Step 5: Commit benchmark results**

```bash
git add tests/benchmark_notebook_performance.sh benchmark_results_*.txt
git commit -m "test: add notebook performance benchmarking

- Script to benchmark Quarto rendering time
- Baseline and optimized results
- Shows speedup from preprocessing optimization"
```

---

## Verification Checklist

Before considering this plan complete, verify:

- [ ] All tests pass: `pytest tests/test_process_variant_table.py -v`
- [ ] Pipeline runs successfully: `snakemake --cores 12 --use-conda`
- [ ] Processed tables exist for all benchmarks
- [ ] Quarto notebook renders without errors
- [ ] Notebook output is identical to original (check figures, tables)
- [ ] Memory usage is reduced (monitor with `htop` during rendering)
- [ ] Rendering time is improved (benchmark results)
- [ ] All commits follow conventional commit format
- [ ] Documentation is updated

## Notes

- This plan uses DRY principle by moving duplicate R logic to reusable Python
- YAGNI: We only process columns actually used in downstream analysis
- TDD: Each task starts with tests before implementation
- Frequent commits after each task completion
- The optimization targets the specific bottleneck (R parsing large tables in memory)
