"""
Rules for generating exclusion intersection tables.

Computes overlap metrics between exclusion BED regions and benchmark regions (dip.bed).
For each benchmark set, generates a table with columns:
- exclusions: Name of the exclusion region
- exclusion_bp: Total bases in the exclusion BED
- intersect_bp: Bases overlapping between exclusion and benchmark regions
- pct_of_dip: Percent of benchmark regions covered by exclusion
"""


# ============================================================================
# Helper Functions
# ============================================================================


def get_exclusion_config(benchmark):
    """Get exclusion configuration for a benchmark set."""
    return config["benchmarksets"].get(benchmark, {}).get("exclusions", [])


def get_exclusion_items(wildcards):
    """Get list of exclusion names for a benchmark set."""
    exclusions = get_exclusion_config(wildcards.benchmark)
    return [item["name"] for item in exclusions]


def get_exclusion_entry(benchmark, exclusion_name):
    """Get a specific exclusion entry by name."""
    exclusions = get_exclusion_config(benchmark)
    for item in exclusions:
        if item["name"] == exclusion_name:
            return item
    raise ValueError(f"Exclusion {exclusion_name} not found for {benchmark}")


def get_exclusion_inputs(wildcards):
    """Get input file paths for an exclusion."""
    entry = get_exclusion_entry(wildcards.benchmark, wildcards.exclusion)
    return [f["path"] for f in entry["files"]]


def get_exclusion_type(wildcards):
    """Get type (single/pair) for an exclusion."""
    entry = get_exclusion_entry(wildcards.benchmark, wildcards.exclusion)
    return entry["type"]


def get_dip_bed_for_exclusions(wildcards):
    """Get the dip.bed path for a benchmark set."""
    return config["benchmarksets"][wildcards.benchmark]["dip_bed"]["path"]


def get_input_checksums(wildcards):
    """Get checksums for exclusion input files."""
    entry = get_exclusion_entry(wildcards.benchmark, wildcards.exclusion)
    return [f["sha256"] for f in entry["files"]]


# ============================================================================
# Rules
# ============================================================================


rule materialize_exclusion:
    """
    Materialize exclusion BED file from source.
    
    For single files: copies the source BED.
    For start/end pairs: concatenates, sorts, and merges into a single BED.
    
    Input checksums are validated in the shell command.
    """
    input:
        files=get_exclusion_inputs,
    output:
        bed=ensure(
            "results/exclusions/{benchmark}/exclusions/{exclusion}.bed",
            non_empty=True,
        ),
    params:
        exclusion_type=get_exclusion_type,
        checksums=get_input_checksums,
    log:
        "logs/exclusions/{benchmark}/{exclusion}_materialize.log",
    message:
        "Materializing exclusion {wildcards.exclusion} for {wildcards.benchmark}"
    threads: 1
    resources:
        mem_mb=2048,
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        echo "Materializing {wildcards.exclusion} ({params.exclusion_type})" > {log}
        echo "Input files: {input.files}" >> {log}
        echo "Started at $(date)" >> {log}
        
        # Validate input checksums
        CHECKSUMS=({params.checksums})
        FILES=({input.files})
        
        for i in "${{!FILES[@]}}"; do
            EXPECTED="${{CHECKSUMS[$i]}}"
            ACTUAL=$(shasum -a 256 "${{FILES[$i]}}" | cut -d ' ' -f 1)
            
            if [ "$EXPECTED" != "$ACTUAL" ]; then
                echo "ERROR: Checksum mismatch for ${{FILES[$i]}}" >> {log}
                echo "  Expected: $EXPECTED" >> {log}
                echo "  Actual:   $ACTUAL" >> {log}
                exit 1
            fi
            echo "Checksum OK: ${{FILES[$i]}}" >> {log}
        done
        
        # Materialize the BED file (always sort for bedtools compatibility)
        if [ "{params.exclusion_type}" == "single" ]; then
            bedtools sort -i {input.files} > {output.bed} 2>> {log}
        else
            # Pair: concatenate, sort, and merge
            cat {input.files} | bedtools sort -i - | bedtools merge -i - > {output.bed} 2>> {log}
        fi
        
        echo "Output lines: $(wc -l < {output.bed})" >> {log}
        echo "Completed at $(date)" >> {log}
        """


rule compute_dip_size:
    """
    Compute total size of benchmark regions (dip.bed) for a benchmark set.
    
    This is computed once per benchmark set and reused for all exclusion
    metrics calculations.
    """
    input:
        bed=get_dip_bed_for_exclusions,
    output:
        size="results/exclusions/{benchmark}/dip_size.txt",
    log:
        "logs/exclusions/{benchmark}/dip_size.log",
    message:
        "Computing dip.bed size for {wildcards.benchmark}"
    threads: 1
    resources:
        mem_mb=2048,
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        echo "Computing dip.bed size for {wildcards.benchmark}" > {log}
        echo "Input: {input.bed}" >> {log}
        echo "Started at $(date)" >> {log}
        
        # Calculate total size (sort then merge overlapping regions)
        bedtools sort -i {input.bed} | bedtools merge -i - | awk '{{sum+=$3-$2}} END {{print sum+0}}' > {output.size} 2>> {log}
        
        echo "Dip size: $(cat {output.size}) bp" >> {log}
        echo "Completed at $(date)" >> {log}
        """


rule compute_exclusion_metrics:
    """
    Compute intersection metrics for an exclusion against the benchmark regions.
    
    Calculates:
    - exclusion_bp: Total bases in the exclusion BED (after merge)
    - intersect_bp: Bases overlapping between exclusion and benchmark regions
    - pct_of_dip: Percent of benchmark regions covered by exclusion
    """
    input:
        exclusion="results/exclusions/{benchmark}/exclusions/{exclusion}.bed",
        dip_bed=get_dip_bed_for_exclusions,
        dip_size="results/exclusions/{benchmark}/dip_size.txt",
    output:
        tsv="results/exclusions/{benchmark}/metrics/{exclusion}.tsv",
    log:
        "logs/exclusions/{benchmark}/{exclusion}_metrics.log",
    message:
        "Computing metrics for {wildcards.exclusion} in {wildcards.benchmark}"
    threads: 1
    resources:
        mem_mb=4096,
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        echo "Computing metrics for {wildcards.exclusion}" > {log}
        echo "Started at $(date)" >> {log}
        
        DIP_SIZE=$(cat {input.dip_size})
        echo "Dip size: $DIP_SIZE bp" >> {log}
        
        # Exclusion size (sort then merge to handle any overlaps)
        EXCL_SIZE=$(bedtools sort -i {input.exclusion} | bedtools merge -i - | awk '{{sum+=$3-$2}} END {{print sum+0}}')
        echo "Exclusion size: $EXCL_SIZE bp" >> {log}
        
        # Intersect size (overlap between exclusion and dip) - sort both inputs
        INTERSECT_SIZE=$(bedtools intersect -a <(bedtools sort -i {input.dip_bed}) -b <(bedtools sort -i {input.exclusion}) | awk '{{sum+=$3-$2}} END {{print sum+0}}')
        echo "Intersect size: $INTERSECT_SIZE bp" >> {log}
        
        # Calculate percentage
        if [ "$DIP_SIZE" -gt 0 ]; then
            PCT=$(awk -v intersect="$INTERSECT_SIZE" -v dip="$DIP_SIZE" 'BEGIN {{printf "%.6f", (intersect / dip) * 100}}')
        else
            PCT=0
        fi
        echo "Percent of dip: $PCT" >> {log}
        
        # Output tab-separated values (no header - added during aggregation)
        echo -e "{wildcards.exclusion}\\t$EXCL_SIZE\\t$INTERSECT_SIZE\\t$PCT" > {output.tsv}
        
        echo "Completed at $(date)" >> {log}
        """


rule aggregate_exclusion_table:
    """
    Aggregate all exclusion metrics into a single CSV table.
    
    Output format:
    - exclusions: Name of the exclusion
    - exclusion_bp: Total bases in exclusion
    - intersect_bp: Overlap with benchmark regions
    - pct_of_dip: Percent of benchmark covered
    """
    input:
        lambda wc: expand(
            "results/exclusions/{benchmark}/metrics/{exclusion}.tsv",
            benchmark=wc.benchmark,
            exclusion=get_exclusion_items(wc),
        ),
    output:
        csv=ensure(
            "results/exclusions/{benchmark}/exclusions_intersection_table.csv",
            non_empty=True,
        ),
    log:
        "logs/exclusions/{benchmark}/aggregate.log",
    message:
        "Aggregating exclusion table for {wildcards.benchmark}"
    threads: 1
    resources:
        mem_mb=1024,
    shell:
        """
        echo "Aggregating exclusion metrics for {wildcards.benchmark}" > {log}
        echo "Input files: {input}" >> {log}
        echo "Started at $(date)" >> {log}
        
        # Write header
        echo "exclusions,exclusion_bp,intersect_bp,pct_of_dip" > {output.csv}
        
        # Concatenate all metric files, converting tabs to commas
        cat {input} | tr '\\t' ',' >> {output.csv}
        
        echo "Output lines: $(wc -l < {output.csv})" >> {log}
        echo "Completed at $(date)" >> {log}
        """
