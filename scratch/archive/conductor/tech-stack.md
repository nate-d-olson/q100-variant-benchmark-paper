# Technology Stack

## Workflow & Pipeline

- **Snakemake:** Primary workflow management system (v9.x) for orchestrating bioinformatics rules and dependencies.
- **Micromamba/Conda:** Used for isolated software environment management via `environment.yaml` and rule-specific environment files in `workflow/envs/`.

## Programming Languages

- **Python (3.14):** Used for custom pipeline scripts, Snakemake logic, and testing (pytest).
- **R (4.5):** Primary language for data analysis, statistical modeling, and plotting (Tidyverse).
- **Bash/Zsh:** Used for shell scripting and automation via the `Makefile`.

## Data Analysis & Reporting

- **Quarto:** Used for generating the reproducible scientific manuscript and supplementary notebooks.
- **Apache Arrow (Parquet):** Core caching and storage format for high-performance data loading in R.
- **Tidyverse:** Collection of R packages for data manipulation and visualization.

## Bioinformatics Core

- **bcftools:** VCF/BCF manipulation and variant annotation.
- **bedtools:** Genome arithmetic and interval analysis.
- **Truvari:** Structural variant comparison and benchmarking.
- **rtg-tools:** Specifically `rtg vcfeval` for small variant comparison.
- **samtools:** SAM/BAM/CRAM processing.

## Quality Assurance & Tooling

- **Linting/Formatting:** `ruff` & `black` (Python), `snakefmt` (Snakemake), `lintr` & `styler` (R), `shellcheck` & `shfmt` (Bash).
- **Testing:** `pytest` (Python/Snakemake scripts), `testthat` (R loading/caching functions).
