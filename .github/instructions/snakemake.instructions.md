---
applyTo: "workflow/**/*.smk,workflow/Snakefile"
---

# Snakemake Rule Writing Standards

When writing or modifying Snakemake rules, follow these mandatory standards:

## Required Directives

Every rule MUST include:

1. **log directive**: Capture stdout/stderr for debugging
   ```python
   log: "logs/{rule_name}/{wildcards}.log"
   ```

2. **ensure() validation**: Validate output file properties
   ```python
   output:
       tsv=ensure("results/{sample}.tsv", non_empty=True)
   ```
   - Use `non_empty=True` for all output files
   - Add `sha256="..."` for deterministic outputs
   - Use `size_gt=1000` for minimum file size requirements

3. **threads directive**: Specify computational resources
   ```python
   threads: 2
   ```

4. **resources directive**: Define memory and other resources
   ```python
   resources:
       mem_mb=4096
   ```

5. **conda directive**: Specify conda environment
   ```python
   conda: "../envs/bcftools.yaml"
   ```

6. **container directive** (optional but recommended):
   ```python
   container: "docker://quay.io/biocontainers/bcftools:1.19"
   ```

## Download Rules

For rules downloading external files:

1. Use `retries: 3` for network stability
2. Wrap outputs in `ancient()` to prevent unnecessary re-downloads
3. Validate checksums from config
4. Log all download activities

Example:
```python
rule download_file:
    output:
        file=ancient("resources/downloads/{dataset}.vcf.gz")
    params:
        url=lambda w: config["downloads"][w.dataset]["url"],
        sha256=lambda w: config["downloads"][w.dataset]["sha256"]
    log:
        "logs/downloads/{dataset}.log"
    retries: 3
    shell:
        """
        wget -O {output.file} {params.url} 2>&1 | tee {log}
        echo "{params.sha256}  {output.file}" | sha256sum -c - 2>&1 | tee -a {log}
        """
```

## Conda Environments

All `workflow/envs/*.yaml` files MUST:
- Use `channels: [conda-forge, bioconda]` (NO `defaults` channel)
- Set `channel_priority: strict`
- Pin specific tool versions (e.g., `bcftools=1.19`)

Example:
```yaml
channels:
  - conda-forge
  - bioconda
channel_priority: strict
dependencies:
  - bcftools=1.19
```

## Wrapper Usage

When using Snakemake wrappers:
- Always use versioned `wrapper_prefix` (set globally in Snakefile)
- Reference wrappers as: `wrapper: f"{wrapper_prefix}/bio/bcftools/query"`
- Check wrapper documentation: https://snakemake-wrappers.readthedocs.io/

## Path Handling

- Use relative paths from rule file: `"../envs/tool.yaml"`
- For Python-based path logic, use `pathlib`:
  ```python
  from pathlib import Path
  output_dir = Path("results") / wildcards.sample
  ```

## Shell Commands

- Use triple-quoted strings for multi-line commands
- Redirect stderr to log: `command 2> {log}` or `2>&1 | tee {log}`
- Escape special characters: `\\t` for tab, `\\n` for newline
- Quote wildcards/params: `{input.vcf}` not `{input.vcf }`

## Best Practices

1. **Descriptive rule names**: Use `generate_sv_len` not `rule1`
2. **Informative log messages**: Include rule context in shell commands
3. **Fail fast**: Exit with error code if validation fails
4. **Comment complex logic**: Explain non-obvious parameter choices
5. **Test locally**: Run `snakemake -n` before committing

## Rule Template

Copy this template when creating new rules. Fill in placeholders and choose the appropriate execution method (shell/script/wrapper):

```python
rule rule_name:
    """
    Brief description of what this rule does.
    
    Detailed explanation of the processing logic, tools used,
    and any important considerations.
    """
    input:
        # Primary input files
        file1="path/to/{wildcard}_input.ext",
        file2="path/to/{wildcard}_input2.ext"
    output:
        # Use ensure() for validation
        result=ensure("results/{wildcard}_output.ext", non_empty=True)
    params:
        # Rule-specific parameters
        param1="value",
        param2=lambda wildcards: config["key"][wildcards.wildcard]
    log:
        "logs/module_name/{wildcard}_rule_name.log"
    threads: 1  # Number of CPU threads
    resources:
        mem_mb=4096,  # Memory in MB
        tmpdir="/tmp"  # Temporary directory
    conda:
        "../envs/tool_name.yaml"
    container:
        "docker://repository/image:tag"  # Optional
    shell:
        """
        # Log execution start
        echo "Processing {wildcards.wildcard}" > {log}
        
        # Main command
        tool_command --option {params.param1} \
            --input {input.file1} \
            --output {output.result} \
            2>> {log}
        
        # Validate output
        if [ ! -s {output.result} ]; then
            echo "ERROR: Output validation failed" >> {log}
            exit 1
        fi
        """
```

### Execution Method Variants

**Choose ONE of these execution methods:**

1. **Shell commands** (default, shown above):
   ```python
   shell:
       """
       command {input} > {output} 2> {log}
       """
   ```

2. **Python/R script**:
   ```python
   script:
       "../scripts/script_name.py"  # or .R, .jl, etc.
   ```
   - Script receives: `snakemake.input`, `snakemake.output`, `snakemake.params`, `snakemake.log`, `snakemake.threads`, `snakemake.resources`, `snakemake.wildcards`

3. **Snakemake wrapper**:
   ```python
   wrapper:
       f"{wrapper_prefix}/bio/tool/command"  # Use global wrapper_prefix
   ```
   - Browse available wrappers: https://snakemake-wrappers.readthedocs.io/
   - Example: `f"{wrapper_prefix}/bio/bcftools/stats"`

### Common Customizations

**Multiple outputs:**
```python
output:
    result1=ensure("results/{sample}_output1.tsv", non_empty=True),
    result2=ensure("results/{sample}_output2.png", size_gt=1000),
    result3=ensure("results/{sample}_output3.vcf.gz", sha256="abc123...")
```

**Conditional parameters:**
```python
params:
    extra=lambda wildcards: "--flag" if wildcards.ref == "GRCh38" else "",
    threshold=lambda wildcards: config["thresholds"][wildcards.vartype]
```

**Dynamic resources:**
```python
resources:
    mem_mb=lambda wildcards, attempt: 4096 * attempt,  # Increase on retry
    runtime=lambda wildcards: 60 if wildcards.ref == "chr21" else 240
```

**Advanced logging:**
```python
shell:
    """
    exec 2> {log}  # Redirect all stderr to log
    set -euo pipefail  # Bash strict mode
    
    echo "Started: $(date)" | tee -a {log}
    command {input} > {output}
    echo "Completed: $(date)" | tee -a {log}
    """
```

## Linting

Before committing, always run:
```bash
snakemake --lint
snakefmt workflow/
```
