# Product Guidelines

## Technical Communication
- **Developer-Centric Style:** Documentation and manuscript drafts should prioritize technical clarity, focusing on the Snakemake rule logic, implementation details, and reproducible code usage.
- **Precision:** Use standard bioinformatics terminology (e.g., "genomic context," "stratification," "liftover") consistently across code, logs, and prose.

## Visual Identity & Figures
- **Publication-Ready Standards:** All figures and tables must meet high-resolution requirements with standardized, colorblind-friendly palettes and consistent typography for academic publication.
- **Traceability:** Maintain a direct, transparent link between the raw pipeline outputs, the R plotting scripts, and the final rendered manuscript artifacts.

## Development & Architecture
- **Snakemake Best Practices:** Adhere strictly to the standard Snakemake distribution structure. Rules must be modular, and environmental dependencies must be explicitly defined via Conda/Mamba.
- **Modular Logic:** Complex processing should reside in dedicated scripts (Python/R) rather than large `run` blocks within the Snakefile.

## Data Integrity & Loading
- **Efficient Caching:** Utilize the project's Parquet-based caching system for all heavy datasets to ensure rapid, reproducible analysis in Quarto notebooks.
- **Schema Validation:** Enforce data consistency by validating all pipeline outputs against the centralized R schema registry.

## Brand & Messaging
- **The Q100 Value Proposition:** Frame the Q100 benchmark as a significant, high-fidelity advancement over the GIAB v4.2.1 release, with a specific emphasis on resolving historically difficult genomic regions.
