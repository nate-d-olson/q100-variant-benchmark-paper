"""
Output organization rules for user-friendly access to final deliverables.

This module creates convenient symlinks to key pipeline outputs, organizing them
by benchmark set rather than by analysis type. This makes it easier for users
to find all outputs related to a specific benchmark.

Original locations remain unchanged - these are just convenient aliases.
"""


rule organize_exclusion_table:
    """
    Create user-friendly symlink to exclusion impact table.

    Provides easy access to exclusion analysis results directly in the
    benchmark directory rather than nested in results/exclusions/.
    """
    input:
        "results/exclusions/{benchmark}/exclusion_impact.csv",
    output:
        "results/benchmarksets/{benchmark}/exclusions.csv",
    log:
        "logs/organize_exclusion_table/{benchmark}.log",
    shell:
        """
        echo "Creating symlink for exclusions table" > {log}
        echo "Source: {input}" >> {log}
        echo "Target: {output}" >> {log}
        ln -sf $(realpath {input}) {output} 2>> {log}
        echo "Symlink created successfully" >> {log}
        """


rule organize_strat_table:
    """
    Create user-friendly symlink to stratification coverage table.

    Provides easy access to stratification metrics directly in the
    benchmark directory rather than nested in results/strat_metrics/.
    """
    input:
        "results/strat_metrics/{benchmark}/stratification_coverage_table.csv",
    output:
        "results/benchmarksets/{benchmark}/stratifications.csv",
    log:
        "logs/organize_strat_table/{benchmark}.log",
    shell:
        """
        echo "Creating symlink for stratifications table" > {log}
        echo "Source: {input}" >> {log}
        echo "Target: {output}" >> {log}
        ln -sf $(realpath {input}) {output} 2>> {log}
        echo "Symlink created successfully" >> {log}
        """


rule organize_variant_table:
    """
    Create user-friendly symlink to variant table.

    Provides easy access to the full variant annotation table directly in the
    benchmark directory rather than nested in results/variant_tables/.
    """
    input:
        "results/variant_tables/{benchmark}/variants.parquet",
    output:
        "results/benchmarksets/{benchmark}/variants.parquet",
    log:
        "logs/organize_variant_table/{benchmark}.log",
    shell:
        """
        echo "Creating symlink for variant table" > {log}
        echo "Source: {input}" >> {log}
        echo "Target: {output}" >> {log}
        ln -sf $(realpath {input}) {output} 2>> {log}
        echo "Symlink created successfully" >> {log}
        """


rule organize_benchmark_outputs:
    """
    Aggregate rule to create all output aliases for a benchmark.

    This convenience rule ensures all user-friendly symlinks are created
    for a given benchmark set. Users can request this rule to get organized
    access to all key outputs.
    """
    input:
        exclusions="results/benchmarksets/{benchmark}/exclusions.csv",
        stratifications="results/benchmarksets/{benchmark}/stratifications.csv",
        variants="results/benchmarksets/{benchmark}/variants.parquet",
    output:
        touch("results/benchmarksets/{benchmark}/.organized"),
    log:
        "logs/organize_benchmark_outputs/{benchmark}.log",
    shell:
        """
        echo "All output aliases created for {wildcards.benchmark}" > {log}
        echo "Available files:" >> {log}
        ls -lh results/benchmarksets/{wildcards.benchmark}/ >> {log} 2>&1
        """
