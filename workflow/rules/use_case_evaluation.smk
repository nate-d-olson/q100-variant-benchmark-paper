"""
SV use-case evaluation.

Wraps `scripts/extract_sv_metrics.py` to aggregate post-refine truvari outputs
from external SV callsets (HiFi/ONT × Sniffles1/Sniffles2) into the three CSVs
consumed by `analysis/use_case_evaluation.qmd`.

Inputs (delivered manually to data/, not produced by this pipeline):
- data/GRCh38_HG002_T~T2TQ100v1.1_Q~{HiFi,ONT}-Snifpub-{sniffles1,sniffles2}_TR~none_GRCh38_HG002-T2TQ100v1.1-dipz2k_stvar-excluded/
  Each contains a truvari bench+refine output subdirectory with
  refine.variant_summary.json. The script discovers the subdir at runtime and
  fails informatively if the parent dir is missing.

Outputs (consumed by analysis/use_case_evaluation.qmd):
- results/use_case/stvar/stratified_metrics.csv
- results/use_case/stvar/svtype_metrics.csv
- results/use_case/stvar/svtype_size_counts.csv

Convenience target: `snakemake use_case_evaluation`
"""

USE_CASE_STVAR_CONTEXTS = ["HP", "TR", "SD", "MAP"]


rule extract_sv_use_case_metrics:
    """Aggregate post-refine truvari outputs into 3 CSVs."""
    input:
        # Stratification BEDs are real pipeline outputs (downloaded by
        # download_stratification). Callset dirs in data/ are externally
        # delivered; the script's existence check handles missing inputs.
        strat_beds=expand(
            "resources/stratifications/GRCh38_{context}.bed.gz",
            context=USE_CASE_STVAR_CONTEXTS,
        ),
        script="scripts/extract_sv_metrics.py",
    output:
        stratified="results/use_case/stvar/stratified_metrics.csv",
        svtype="results/use_case/stvar/svtype_metrics.csv",
        size_counts="results/use_case/stvar/svtype_size_counts.csv",
    log:
        "logs/extract_sv_use_case_metrics.log",
    conda:
        "../envs/truvari.yaml"
    shell:
        "python {input.script} > {log} 2>&1"


rule use_case_evaluation:
    """Convenience target: `snakemake use_case_evaluation`."""
    input:
        rules.extract_sv_use_case_metrics.output,
