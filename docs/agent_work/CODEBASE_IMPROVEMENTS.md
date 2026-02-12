# Q100 Variant Benchmark Pipeline - Suggested Improvements

This document outlines 5 key improvements identified for the q100-variant-benchmark codebase to enhance maintainability, reliability, and usability.

---

## 1. Fix Critical Bug in `stratify_comparison.py`

**Priority: High**
**Impact: Critical - Script will fail**

### Problem

Line 125 in `workflow/scripts/stratify_comparison.py` has incorrect indentation - `import gzip` is placed inside the `main()` function instead of at the module level. This breaks the script's functionality.

```python
# Current (INCORRECT):
def main():
    logging.info("Starting stratification analysis")

import gzip    # <-- This should be at the top with other imports
```

### Solution

Move `import gzip` to the top of the file with other imports (lines 1-5).

### Benefits

- Fixes immediate script execution error
- Follows Python best practices for import placement
- Improves code readability

---

## 2. Add Comprehensive Test Suite

**Priority: High**
**Impact: Reliability and maintainability**

### Problem

The pipeline currently has no automated tests for:

- Python helper scripts (`workflow/scripts/*.py`)
- Snakemake rule logic (`workflow/rules/common.smk` functions)
- Configuration validation beyond schema

### Solution

Implement pytest-based testing:

```bash
tests/
├── unit/
│   ├── test_common_helpers.py          # Test helper functions
│   ├── test_stratify_comparison.py     # Test new comparison script
│   └── test_extract_info_fields.py     # Test VCF parsing
├── integration/
│   ├── test_downloads.py               # Test download rules
│   └── test_comparisons_workflow.py    # Test benchmark comparisons
├── fixtures/
│   ├── small_test.vcf.gz               # Minimal test VCF
│   └── test_config.yaml                # Test configuration
└── conftest.py                          # Pytest configuration
```

**Key test areas:**

1. Helper function logic (common.smk)
2. VCF parsing and annotation scripts
3. Stratification analysis correctness
4. Configuration validation
5. Integration tests for critical workflows

### Benefits

- Catch bugs before they reach production
- Enable confident refactoring
- Document expected behavior
- Enable CI/CD integration

---

## 3. Improve Error Handling and Validation

**Priority: Medium**
**Impact: User experience and debugging**

### Problem

Current error handling is minimal:

- `stratify_comparison.py` logs errors but continues processing
- No pre-flight validation of inputs (URLs, file existence)
- Generic subprocess error messages lack context
- Configuration errors surface late in pipeline execution

### Solution

**A. Add Pre-flight Validation Rule:**

```python
rule validate_inputs:
    """Validate all inputs before pipeline execution"""
    input:
        config="config/config.yaml"
    output:
        touch("results/.validation_complete")
    script:
        "scripts/validate_inputs.py"
```

**B. Enhance `stratify_comparison.py` error handling:**

```python
# Instead of logging and returning 0, raise informative exceptions
def run_bedtools_intersect_vcf(vcf_file, bed_file):
    if not os.path.exists(vcf_file):
        raise FileNotFoundError(f"VCF file not found: {vcf_file}")
    if not os.path.exists(bed_file):
        raise FileNotFoundError(f"BED file not found: {bed_file}")

    try:
        result = subprocess.run(
            ["bedtools", "intersect", "-u", "-a", vcf_file, "-b", bed_file],
            capture_output=True,
            text=True,
            check=True  # Raise on non-zero exit
        )
        # ... process output
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"bedtools intersect failed:\n"
            f"  VCF: {vcf_file}\n"
            f"  BED: {bed_file}\n"
            f"  Error: {e.stderr}"
        ) from e
```

**C. Add input validation script:**

- Check URLs are reachable (HTTP HEAD request)
- Verify required tools exist (bedtools, rtg, truvari)
- Validate stratification file existence
- Check disk space for large downloads

### Benefits

- Fail fast with clear error messages
- Reduce debugging time
- Prevent partial pipeline execution
- Better user experience

---

## 4. Enhance Documentation

**Priority: Medium**
**Impact: Usability and maintainability**

### Problem

Current documentation gaps:

- No documentation for new comparison functionality (Gemini addition)
- Missing docstrings in `stratify_comparison.py`
- No type hints in several helper functions
- README doesn't explain comparison outputs
- No examples of how to add new comparisons

### Solution

**A. Add Comparison Documentation Section to README.qmd:**

```markdown
## Benchmark Comparisons

The pipeline supports comparing different benchmark versions to quantify:
- Shared and unique variants between versions
- Shared and unique genomic regions (bases)
- Stratification-specific differences

### Configured Comparisons

| Comparison ID | Description |
|---------------|-------------|
| v5_vs_v4_smvar | v5.0q vs v4.2.1 small variants (GRCh38) |
| v5_vs_v0.6_stvar | v5.0q vs v0.6 structural variants (GRCh37) |

### Outputs

- `results/stats/{comp_id}_variants.csv` - Variant counts (shared/unique)
- `results/stats/{comp_id}_regions.csv` - Base pair counts (shared/unique)

### Adding New Comparisons

Edit `config/config.yaml`:
```yaml
comparisons:
  my_comparison:
    type: smvar  # or stvar
    new_benchmark: v5.0q_GRCh38_smvar
    old_benchmark: v4.2.1_GRCh38_smvar
    ref: GRCh38
```

```

**B. Add comprehensive docstrings to `stratify_comparison.py`:**
```python
def run_bedtools_intersect_vcf(vcf_file: str, bed_file: str) -> int:
    """
    Count variants in VCF overlapping with BED regions.

    Uses bedtools intersect with -u flag to count unique overlapping
    variant records (excluding duplicates).

    Args:
        vcf_file: Path to input VCF file (can be gzipped)
        bed_file: Path to BED file defining regions of interest

    Returns:
        Number of variants overlapping with BED regions

    Raises:
        FileNotFoundError: If input files don't exist
        subprocess.CalledProcessError: If bedtools fails

    Example:
        >>> count = run_bedtools_intersect_vcf("vars.vcf.gz", "regions.bed")
        >>> print(f"Found {count} variants in regions")
    """
```

**C. Add type hints to `common.smk` functions:**

```python
from typing import List, Dict, Any

def get_exclusion_file_path(benchmark: str, exclusion_name: str, file_idx: int) -> str:
    """Get the standardized local path for an exclusion file."""
    return f"resources/exclusions/{benchmark}/{exclusion_name}_{file_idx}.bed"
```

### Benefits

- Easier onboarding for new contributors
- Self-documenting code
- IDE autocomplete support
- Reduced knowledge silos

---

## 5. Optimize `stratify_comparison.py` Performance

**Priority: Low**
**Impact: Pipeline execution speed**

### Problem

Current implementation has performance inefficiencies:

- Multiple subprocess calls per stratification (N × 3 bedtools calls)
- Temporary files created/destroyed repeatedly
- Sequential processing (no parallelization)
- Gzipped files decompressed on every access

### Solution

**A. Batch bedtools operations:**

```python
# Instead of N separate bedtools calls, use bedtools multiIntersect
def batch_intersect_stratifications(vcf, strat_beds):
    """Process all stratifications in a single bedtools call."""
    # Create temporary file with all stratification beds
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt') as bedlist:
        bedlist.write('\n'.join(strat_beds))
        bedlist.flush()

        # Use bedtools multiIntersect for efficient batch processing
        result = subprocess.run(
            ["bedtools", "multiIntersect", "-i", bedlist.name],
            capture_output=True, text=True, check=True
        )
        # Parse multiIntersect output...
```

**B. Cache decompressed files:**

```python
@functools.lru_cache(maxsize=32)
def get_decompressed_vcf(vcf_path: str) -> str:
    """Cache decompressed VCFs for repeated access."""
    if not vcf_path.endswith('.gz'):
        return vcf_path

    cache_path = f"/tmp/{os.path.basename(vcf_path)}.cache"
    if not os.path.exists(cache_path):
        with gzip.open(vcf_path, 'rb') as f_in:
            with open(cache_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    return cache_path
```

**C. Use pybedtools instead of subprocess:**

```python
import pybedtools

def count_overlaps(vcf_file: str, bed_file: str) -> int:
    """Count overlaps using pybedtools (faster than subprocess)."""
    vcf = pybedtools.BedTool(vcf_file)
    regions = pybedtools.BedTool(bed_file)
    return len(vcf.intersect(regions, u=True))
```

### Benefits

- Reduce pipeline runtime (estimated 30-50% for comparison rules)
- Lower disk I/O
- More Pythonic code (easier to maintain)
- Better resource utilization

---

## Implementation Priority

1. **Fix Bug (Improvement #1)** - Immediate
2. **Add Tests (Improvement #2)** - Next sprint
3. **Improve Error Handling (Improvement #3)** - Next sprint
4. **Documentation (Improvement #4)** - Ongoing
5. **Performance Optimization (Improvement #5)** - Future enhancement

---

## Success Metrics

- **Bug Fix:** Script executes without import errors
- **Tests:** >80% code coverage, all tests passing in CI
- **Error Handling:** Zero generic "subprocess failed" errors in logs
- **Documentation:** All public functions have docstrings with type hints
- **Performance:** Comparison rules complete 30%+ faster

---

*Document created: 2026-01-13*
*Last updated: 2026-01-13*
