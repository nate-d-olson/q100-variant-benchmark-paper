---
applyTo: "workflow/scripts/**/*.py"
---

# Python Snakemake Scripts Instructions

This document provides guidelines for Python scripts used in Snakemake rules.

## Script Purpose

Python scripts in `workflow/scripts/` are invoked by Snakemake rules and receive parameters through the special `snakemake` object:
- `snakemake.input` - Input files
- `snakemake.output` - Output files
- `snakemake.params` - Parameters from rule
- `snakemake.log` - Log file path
- `snakemake.threads` - Thread count
- `snakemake.resources` - Resource specifications
- `snakemake.wildcards` - Wildcard values

## Code Style

**Formatter:** `ruff` (configured in Makefile)
- Run formatting: `ruff format workflow/scripts/`
- Check formatting: `ruff check workflow/scripts/`

**Linter:** `ruff`
- Configuration: Special per-file-ignores in pyproject.toml (if exists)
- Snakemake entrypoint scripts: F821 (undefined names) allowed due to `snakemake` object
- Utility modules: F821 NOT allowed (must detect undefined names)

## Script Template

```python
#!/usr/bin/env python3
# AI Disclosure: This script was developed with assistance from [AI Tool Name].
"""
Brief description of script purpose.

Detailed explanation of inputs, outputs, and processing logic.
"""

import sys
from pathlib import Path

def main():
    """Main script logic."""
    # Access Snakemake parameters
    input_file = snakemake.input[0]
    output_file = snakemake.output[0]
    log_file = snakemake.log[0]
    
    try:
        # Open log file for writing
        with open(log_file, 'w') as log:
            log.write(f"Processing {input_file}\\n")
            
            # Script logic here
            process_data(input_file, output_file, log)
            
            log.write("Processing complete\\n")
    
    except Exception as e:
        with open(log_file, 'a') as log:
            log.write(f"ERROR: {str(e)}\\n")
        raise

def process_data(input_file, output_file, log):
    """Process data (implementation)."""
    pass

if __name__ == "__main__":
    main()
```

## Resource Management

**ALWAYS use context managers (with statement) for file handles:**

```python
# Good - ensures file is closed even if error occurs
with open(file_path, 'r') as f:
    data = f.read()

# Bad - file may not be closed on exception
f = open(file_path, 'r')
data = f.read()
f.close()
```

This pattern is consistently used across all Python scripts in this project.

## Error Handling

**Snakemake Python scripts should raise exceptions (not return empty values) when operations fail:**

```python
# Good - ensures workflow fails visibly
result = subprocess.run(cmd, capture_output=True, text=True)
if result.returncode != 0:
    raise RuntimeError(f"Command failed: {result.stderr}")

# Bad - silently returns empty result
result = subprocess.run(cmd, capture_output=True, text=True)
if result.returncode != 0:
    return {}  # Don't do this!
```

This ensures the Snakemake workflow fails fast and visibly when errors occur.

## Common Patterns

### Truvari Integration

When working with variant data, use Truvari's VariantRecord API:

```python
import truvari

# Read VCF with VariantFile
vcf = truvari.VariantFile(vcf_path)

for entry in vcf:
    vr = truvari.VariantRecord(entry)
    
    # Get variant type (returns SV enum)
    var_type = vr.var_type().name  # Use .name to extract string
    
    # Get variant size
    var_size = vr.var_size()
    
    # Get genotype (returns GT enum)
    gt = truvari.get_gt(vr.gt()).name  # Use .name to extract string
    
    # Get size bin (returns string directly)
    szbin = truvari.get_sizebin(var_size)
```

**Important:** Use `.name` to extract strings from Truvari enums:
- `vr.var_type()` returns SV enum → use `.name` to get "DEL", "INS", "SNP"
- `truvari.get_gt()` returns GT enum → use `.name` to get "HET", "HOM", "REF"
- `truvari.get_sizebin()` returns string directly (no enum)
- **NEVER** use `str()` on enums - gives wrong format like "SV.DEL" instead of "DEL"

### CSV/TSV I/O

Use pandas for structured data:

```python
import pandas as pd

# Read CSV/TSV
df = pd.read_csv(input_file, sep='\\t')

# Write CSV/TSV
df.to_csv(output_file, sep='\\t', index=False)
```

### Parquet I/O

Use PyArrow or pandas for Parquet:

```python
import pyarrow.parquet as pq
import pandas as pd

# Write with PyArrow (more control)
table = pa.Table.from_pandas(df, schema=schema)
pq.write_table(table, output_file, compression='zstd')

# Write with pandas (simpler)
df.to_parquet(output_file, compression='zstd', index=False)
```

### BED File Processing

Use bedtools via subprocess or pybedtools:

```python
import subprocess

# Using bedtools via subprocess
cmd = [
    'bedtools', 'intersect',
    '-a', input_a,
    '-b', input_b,
    '-wa', '-wb'
]
result = subprocess.run(cmd, capture_output=True, text=True, check=True)

# Or use pybedtools (if available)
import pybedtools
a = pybedtools.BedTool(input_a)
b = pybedtools.BedTool(input_b)
result = a.intersect(b, wa=True, wb=True)
result.saveas(output_file)
```

### Logging Best Practices

Always log important steps and errors:

```python
with open(snakemake.log[0], 'w') as log:
    log.write(f"Starting processing at {datetime.now()}\\n")
    log.write(f"Input: {input_file}\\n")
    log.write(f"Output: {output_file}\\n")
    
    # Processing steps
    log.write(f"Processing {count} records...\\n")
    
    # Results
    log.write(f"Wrote {output_count} records\\n")
    log.write(f"Completed at {datetime.now()}\\n")
```

## Script Categories

### Data Processing Scripts
- `generate_variant_parquet.py` - Convert VCF to Parquet using Truvari
- `count_variants_by_genomic_context.py` - Count variants per genomic context
- `combine_beds_with_id.py` - Merge BED files with ID tracking

### Analysis Scripts
- `compute_bed_metrics.py` - Calculate BED region statistics
- `compute_exclusion_interactions.py` - Upset-style exclusion analysis
- `annotate_old_benchmark_status.py` - Cross-version comparison

### Utility Scripts
- `exceptions.py` - Custom exception classes
- `generate_header_lines.py` - Generate VCF header annotations

## Dependencies

Common Python packages used:
- **pandas** - Data manipulation
- **numpy** - Numerical operations
- **pyarrow** - Parquet I/O
- **truvari** - Variant analysis utilities
- **pybedtools** - BED file operations (optional)
- **pysam** - SAM/BAM/VCF I/O (via Truvari)

## Testing

Python tests use pytest:

```bash
# Run all Python tests
pytest tests/

# Run specific test file
pytest tests/test_scripts.py

# Run with coverage
pytest --cov=workflow/scripts tests/
```

## Before Committing

1. Format Python code: `ruff format workflow/scripts/`
2. Check linting: `ruff check workflow/scripts/`
3. Run tests: `pytest tests/`
4. Test script in Snakemake workflow: `snakemake -n --quiet`

## Common Issues

### Tuple/String Conversion in Parquet

When VCF INFO fields with multiple values are read by Truvari and written to Parquet, they can become string representations of tuples:
- Problem: `"('HP', 'MAP')"` instead of `['HP', 'MAP']`
- Solution: Convert tuples to lists before writing to Parquet

```python
# Convert tuple-like INFO fields to lists
if isinstance(context_ids, (tuple, list)):
    context_ids = list(context_ids)
```

### Stale Output Cache

Snakemake does not automatically detect when a Python script changes:
- When you modify a script, delete the affected output directories
- This forces Snakemake to regenerate outputs with the updated script

```bash
# Example: after modifying generate_variant_parquet.py
rm -rf results/variant_tables/
snakemake --cores 4
```

## Script Execution

Scripts are executed by Snakemake with the following pattern:

```python
script:
    "../scripts/script_name.py"
```

The script receives the `snakemake` object automatically and should:
1. Read inputs from `snakemake.input`
2. Write outputs to `snakemake.output`
3. Log to `snakemake.log`
4. Use `snakemake.threads` for parallelization
5. Access wildcards via `snakemake.wildcards`
6. Access params via `snakemake.params`
