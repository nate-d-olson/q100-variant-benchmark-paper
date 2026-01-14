# Q100 Variant Benchmark Pipeline - Troubleshooting Guide

## Common Issues and Solutions

### Download Failures

#### Issue: SHA256 checksum mismatch
```
Error: SHA256 checksum verification failed
Expected: abc123...
Actual: def456...
```

**Causes:**
- Incomplete download due to network interruption
- Corrupted transfer
- File changed at source

**Solutions:**
1. Delete the partially downloaded file:
   ```bash
   rm resources/benchmarks/{benchmark}/{ref}/variants.vcf.gz
   ```

2. Re-run the pipeline:
   ```bash
   snakemake --cores 4
   ```

3. If checksum continues to fail, verify the URL and checksum in `config/config.yaml`

---

#### Issue: Connection timeout during download
```
Error: Read timed out after 300 seconds
```

**Causes:**
- Slow or unstable network connection
- Large file download interrupted

**Solutions:**
1. Increase timeout in download rules (edit `workflow/rules/downloads.smk`):
   ```python
   params:
       timeout=3600  # Increase from 300 to 3600 seconds
   ```

2. Use manual download with resume capability:
   ```bash
   wget -c -O resources/benchmarks/v5q/GRCh38/variants.vcf.gz [URL]
   ```

3. Verify downloaded file checksum manually:
   ```bash
   sha256sum resources/benchmarks/v5q/GRCh38/variants.vcf.gz
   ```

---

### Memory Issues

#### Issue: Killed process during VCF annotation
```
Killed
[snakemake] Error in rule annotate_vcf
```

**Causes:**
- Insufficient memory for large VCF files (multi-GB files)
- Too many parallel jobs competing for memory

**Solutions:**
1. Reduce parallel jobs:
   ```bash
   snakemake --cores 1  # Run sequentially
   ```

2. Request more memory (if using cluster):
   ```python
   # Add to rule
   resources:
       mem_mb=32000  # Request 32GB
   ```

3. Monitor memory usage:
   ```bash
   watch -n 1 'free -h'
   ```

---

#### Issue: Large stratification TSV files (600MB+)
```
Warning: v5_grch38_strat.tsv is 650MB
```

**Causes:**
- Comprehensive stratification across entire genome
- Normal for this pipeline

**Solutions:**
1. This is expected behavior - no action needed

2. If disk space is limited, add cleanup rules:
   ```python
   temp("results/var_tables/{benchmark}/{ref}/annotated.tsv")
   ```

3. Use `--delete-temp-output` flag to auto-clean temporary files:
   ```bash
   snakemake --cores 4 --delete-temp-output
   ```

---

### Conda Environment Conflicts

#### Issue: Package version conflicts
```
Error: Solving environment: failed
Conflicts: bcftools-1.18 requires libdeflate>=1.10
```

**Causes:**
- Incompatible package versions in environment files
- Outdated conda solver

**Solutions:**
1. Update conda/mamba:
   ```bash
   mamba update -n base mamba
   ```

2. Clear conda cache:
   ```bash
   mamba clean --all
   ```

3. Force environment recreation:
   ```bash
   snakemake --cores 4 --use-conda --conda-create-envs-only --force
   ```

4. Pin specific versions in `workflow/envs/*.yaml`:
   ```yaml
   dependencies:
     - bcftools=1.17  # Pin specific version
   ```

---

#### Issue: Conda environment activation fails
```
Error: CondaEnvException: Could not find conda environment: /path/to/env
```

**Causes:**
- Corrupted environment
- Missing conda/mamba installation

**Solutions:**
1. Verify conda/mamba is in PATH:
   ```bash
   which mamba
   ```

2. Recreate environments:
   ```bash
   rm -rf .snakemake/conda
   snakemake --cores 4 --use-conda --conda-create-envs-only
   ```

3. Use micromamba if mamba fails:
   ```bash
   snakemake --cores 4 --use-conda --conda-frontend micromamba
   ```

---

### Data Validation Errors (NEW)

#### Issue: VCF header validation fails
```
DataFormatError: VCF header missing ##fileformat line
File: resources/benchmarks/v5q/GRCh38/variants.vcf.gz
```

**Causes:**
- Corrupted VCF file
- Incomplete download
- Non-VCF file with wrong extension

**Solutions:**
1. Check file is actually a VCF:
   ```bash
   zcat resources/benchmarks/v5q/GRCh38/variants.vcf.gz | head -20
   ```

2. Re-download the file (see download troubleshooting above)

3. Validate with bcftools:
   ```bash
   bcftools view --header-only resources/benchmarks/v5q/GRCh38/variants.vcf.gz
   ```

---

#### Issue: BED file not sorted
```
ValidationError: File not sorted - position 10500 comes after 10600 on chr1
Line: 523
```

**Causes:**
- BED file not sorted by position
- Custom BED file added without sorting

**Solutions:**
1. Sort BED file:
   ```bash
   sort -k1,1 -k2,2n input.bed > sorted.bed
   mv sorted.bed input.bed
   ```

2. Use bedtools sort:
   ```bash
   bedtools sort -i input.bed > sorted.bed
   ```

3. Disable sort checking (if intentional):
   - Edit `workflow/scripts/validators.py`:
   ```python
   validate_bed_format(file_path, check_sorted=False)
   ```

---

### Pipeline Execution Errors

#### Issue: Missing input files
```
MissingInputException: Missing input files for rule all:
resources/benchmarks/v5q/GRCh38/variants.vcf.gz
```

**Causes:**
- Download rule hasn't run
- File path mismatch in config
- Download failed silently

**Solutions:**
1. Run downloads explicitly:
   ```bash
   snakemake --cores 4 --forceall download_benchmark_vcf
   ```

2. Check config URLs are accessible:
   ```bash
   curl -I [URL from config.yaml]
   ```

3. Verify file paths match config:
   ```bash
   grep -A 5 "v5q:" config/config.yaml
   ```

---

#### Issue: Locked directory
```
WorkflowError: Directory cannot be locked
```

**Causes:**
- Previous Snakemake run crashed
- Lock file not cleaned up

**Solutions:**
1. Remove lock file:
   ```bash
   rm -rf .snakemake/locks
   ```

2. Force unlock:
   ```bash
   snakemake --unlock
   ```

3. Check for zombie processes:
   ```bash
   ps aux | grep snakemake
   kill [PID if found]
   ```

---

### Logging and Debugging

#### Issue: Need more detailed logs
```
Rule X failed but log doesn't show why
```

**Solutions:**
1. Enable verbose logging:
   ```bash
   snakemake --cores 4 --verbose
   ```

2. Check rule-specific log files:
   ```bash
   cat logs/{rule_name}/{wildcards}.log
   ```

3. Run rule in isolation with detailed output:
   ```bash
   snakemake --cores 1 --forcerun {rule_name} --printshellcmds
   ```

4. For Python scripts with new logging (NEW):
   ```bash
   # Logs now include structured format:
   # [2026-01-13 10:30:45] [INFO] [script_name] Message
   tail -f logs/{rule_name}/{wildcards}.log
   ```

---

#### Issue: Understanding pipeline progress
```
How do I know if my pipeline is stuck or just slow?
```

**Solutions:**
1. Use `--verbose` to see active jobs:
   ```bash
   snakemake --cores 4 --verbose
   ```

2. Monitor system resources:
   ```bash
   # Terminal 1
   htop

   # Terminal 2
   watch -n 2 'ls -lh results/*/* | tail -20'
   ```

3. Generate DAG to understand dependencies:
   ```bash
   snakemake --dag | dot -Tpdf > dag.pdf
   ```

4. Use dry-run to see what will execute:
   ```bash
   snakemake --cores 4 --dry-run
   ```

---

### Configuration Errors

#### Issue: Invalid YAML syntax
```
SyntaxError: while scanning a simple key
  in "config/config.yaml", line 42, column 3
```

**Causes:**
- Incorrect YAML indentation
- Missing quotes around special characters
- Tabs instead of spaces

**Solutions:**
1. Validate YAML syntax:
   ```bash
   python -c "import yaml; yaml.safe_load(open('config/config.yaml'))"
   ```

2. Use YAML linter:
   ```bash
   yamllint config/config.yaml
   ```

3. Common fixes:
   - Use spaces, not tabs
   - Quote strings with special characters: `url: "http://example.com?query=value"`
   - Ensure consistent indentation (2 or 4 spaces)

---

#### Issue: Schema validation fails
```
WorkflowError: Config validation failed
```

**Causes:**
- Missing required fields
- Incorrect data types
- Invalid values

**Solutions:**
1. Check schema requirements:
   ```bash
   cat config/schema/config.schema.yaml
   ```

2. Validate config against schema manually:
   ```bash
   python -c "from snakemake.utils import validate; import yaml; config=yaml.safe_load(open('config/config.yaml')); validate(config, 'config/schema/config.schema.yaml')"
   ```

3. Compare with example config in documentation

---

## Performance Optimization

### Slow Pipeline Execution

**Issue**: Pipeline takes hours to complete

**Diagnostic Steps**:
1. Identify bottleneck rules:
   ```bash
   snakemake --cores 4 --profile --printshellcmds 2>&1 | grep "seconds"
   ```

2. Check if rules can be parallelized:
   ```bash
   snakemake --dag | grep -A 5 -B 5 "bottleneck_rule"
   ```

**Solutions**:
1. Increase parallel jobs:
   ```bash
   snakemake --cores 8  # Use more cores
   ```

2. Use cluster execution for large datasets:
   ```bash
   snakemake --cluster "sbatch --mem=16G --time=4:00:00" --jobs 20
   ```

3. Add benchmark directives to rules (see `docs/performance.md` - coming soon)

---

### Disk Space Issues

**Issue**: Running out of disk space

**Diagnostic**:
```bash
du -sh results/* | sort -h
```

**Solutions**:
1. Clean intermediate files:
   ```bash
   snakemake --delete-temp-output
   ```

2. Remove old logs:
   ```bash
   find logs -mtime +30 -delete  # Delete logs older than 30 days
   ```

3. Compress large TSV files:
   ```bash
   gzip results/var_tables/*/*.tsv
   ```

---

## Getting Help

### Reporting Issues

When reporting issues, include:

1. **Exact error message**:
   ```bash
   snakemake --cores 4 2>&1 | tee error.log
   ```

2. **Snakemake version**:
   ```bash
   snakemake --version
   ```

3. **Conda environment info**:
   ```bash
   conda list
   ```

4. **Config file** (sanitize sensitive data):
   ```bash
   cat config/config.yaml
   ```

5. **Log files** from failed rules:
   ```bash
   cat logs/{rule_name}/{wildcards}.log
   ```

### Useful Diagnostic Commands

```bash
# Show pipeline status
snakemake --summary

# List all rules
snakemake --list

# Show reason for rule execution
snakemake --reason

# Print shell commands without executing
snakemake --dry-run --printshellcmds

# Generate rule graph
snakemake --rulegraph | dot -Tpng > rulegraph.png

# Check for circular dependencies
snakemake --dag 2>&1 | grep -i "cycle"
```

---

## Quick Reference

| Symptom | Likely Cause | Quick Fix |
|---------|-------------|-----------|
| "Checksum mismatch" | Incomplete download | Delete file, re-run |
| "Killed" | Out of memory | Reduce `--cores` |
| "Missing input" | Download failed | Check logs, re-run downloads |
| "Locked directory" | Crashed run | `snakemake --unlock` |
| "Invalid YAML" | Syntax error | Check indentation, quotes |
| Slow execution | Single-threaded | Increase `--cores` |
| Disk full | Temp files | `--delete-temp-output` |
| Conda error | Corrupted env | `rm -rf .snakemake/conda` |

---

*Last Updated: 2026-01-13*
*Branch: feature/codebase-improvements*
