.PHONY: help dry-run lint format format-check test test-py test-r clean clean-deep dag run pre-commit-install chr8-preview ideogram ideogram-karyoscope

QMD_FILES := $(shell find . -name '*.qmd' -not -path './.snakemake/*' -not -path './results/*' -not -path './logs/*')
MD_FILES := "**/*.md"

# Default target
help:
	@echo "Q100 Variant Benchmark Analysis - Makefile"
	@echo ""
	@echo "Available targets:"
	@echo "  help             - Show this help message"
	@echo "  dry-run          - Validate workflow with dry-run"
	@echo "  lint             - Run linting (Snakemake + Python + R/Quarto)"
	@echo "  format           - Format all code (Python, Snakemake, Markdown, R)"
	@echo "  format-check     - Check formatting without modifying files"
	@echo "  test             - Run all tests (lint + format-check + dry-run)"
	@echo "  test-py          - Run self-contained Python unit tests"
	@echo "  test-r           - Run R unit tests"
	@echo "  pre-commit-install - Install pre-commit hooks for automatic formatting"
	@echo "  dag              - Generate pipeline DAG visualization (PDF + DOT)"
	@echo "  run              - Execute the pipeline with conda environments"
	@echo "  clean            - Remove logs and temporary files"
	@echo "  clean-deep       - Remove heavier local build/test/cache artifacts"
	@echo "  chr8-preview     - Regenerate chr8 figure and open PDF at actual size"
	@echo ""
	@echo "Usage:"
	@echo "  make format          # format all files"
	@echo "  make lint            # lint all files"
	@echo "  make test            # run all checks"
	@echo "  make pre-commit-install  # set up git hooks"

# Dry-run workflow validation
dry-run:
	@echo "==> Running workflow dry-run validation..."
	snakemake -n --quiet

# Lint workflow
lint-smk:
	@echo "==> Linting Snakemake workflow..."
	snakemake --lint
lint-py:
	@echo "==> Linting Python scripts..."
	ruff check workflow/scripts/
lint-r:
	@echo "==> Linting R scripts and Quarto files..."
	Rscript -e "files <- list.files('.', pattern = '.R$$', full.names = TRUE); paths <- c('R', 'tests','analysis'); files <- c(files, unlist(lapply(paths[dir.exists(paths)], function(p) list.files(p, pattern = '(R|r|qmd)$$', recursive = TRUE, full.names = TRUE)))); lints <- do.call(c, lapply(files, lintr::lint)); if (length(lints) > 0) { print(lints); quit(status = 1) }"
lint-md:
	@echo "==> Linting Markdown files..."
	markdownlint $(MD_FILES)

lint: lint-smk lint-py lint-r

# Format all code files
format:
	@echo "==> Formatting Python scripts..."
	ruff format workflow/scripts/
	ruff check --select=F --ignore=F821 --fix workflow/scripts/
	@echo "==> Formatting Snakemake files..."
	snakefmt workflow/
	@echo "==> Formatting Markdown files..."
	pre-commit run prettier --all-files || true
	@echo "==> Formatting R files..."
	Rscript -e 'paths <- c("R","scripts","analysis"); paths <- paths[dir.exists(paths)]; if (length(paths) > 0) { lapply(paths, function(p) styler::style_dir(p, recursive = TRUE)) }' || echo "(R/styler not available, skipping)"

format-check:
	@echo "==> Checking Python formatting..."
	ruff format --check workflow/scripts/
	@echo "==> Checking Snakefile formatting..."
	snakefmt --check workflow/
	@echo "==> Checking Markdown formatting..."
	pre-commit run prettier --all-files
	@echo "==> Checking R formatting..."
	Rscript -e 'paths <- c("R","scripts"); paths <- paths[dir.exists(paths)]; if (length(paths) > 0) { lapply(paths, function(p) styler::style_dir(p, recursive = TRUE, dry = "fail")) }' || echo "(R/styler not available, skipping)"

# Install pre-commit hooks
pre-commit-install:
	pip install pre-commit
	pre-commit install
	@echo "==> Pre-commit hooks installed. Formatting will run automatically on git commit."

# Run self-contained unit tests. test_common_helpers.py targets removed
# Snakemake helpers and is excluded until it is rewritten or removed.
test-py:
	pytest tests/unit/ --ignore=tests/unit/test_common_helpers.py -v

test-r:
	Rscript -e 'files <- list.files("tests", pattern = "^test_.*[.]R$$", full.names = TRUE); stopifnot(length(files) > 0L); for (file in files) { message("Running ", file); source(file) }'

# Run all checks
test: lint format-check test-py test-r dry-run
	@echo "==> All tests passed!"

# Generate pipeline DAG visualization
dag:
	@echo "==> Generating pipeline DAG..."
	@mkdir -p results/dag
	snakemake --dag > results/dag/pipeline.dot
	@echo "==> Machine-readable DOT: results/dag/pipeline.dot"
	@echo "==> Generating PDF visualization..."
	dot -Tpdf results/dag/pipeline.dot > results/dag/pipeline.pdf
	echo "==> PDF visualization: results/dag/pipeline.pdf"
	@echo "==> DAG generation complete"

# Run the pipeline
run:
	@echo "==> Running pipeline with conda environments..."
	time snakemake --cores 20 --sdm conda --conda-frontend conda --report pipeline_run.html --report-after-run
	@echo "==> Pipeline execution complete"

# Regenerate chr8 synteny figure and open the PDF at actual manuscript size.
# PDFs carry physical dimensions — Preview → View → Actual Size shows true inches.
# Edit CONFIG constants in workflow/scripts/make_chr8_figure.py, then rerun.
chr8-preview:
	@echo "==> Regenerating chr8 synteny figure..."
	rm -f results/chr8_synteny/chr8_figure.pdf results/chr8_synteny/chr8_figure.png \
	      results/chr8_synteny/plotsr_chr8_full.pdf results/chr8_synteny/plotsr_chr8_zoom.pdf
	snakemake --sdm conda --cores 4 chr8_synteny
	@echo "==> Opening PDF (View → Actual Size for true manuscript dimensions)..."
	open results/chr8_synteny/chr8_figure.pdf

# Generate genome-wide ideogram figure
ideogram: scripts/make_ideogram.R resources/hg19ToHg38.over.chain.gz
	Rscript scripts/make_ideogram.R

# Alternative genome-view ideogram rendered with KaryoScope's painted-chromosome
# renderer: horizontal (landscape) chromosomes painted with the actual benchmark
# regions, with telomere markers (figures/ideogram_karyoscope.{svg,pdf,png}).
# Requires the karyoscope conda env (one-time):
#   mamba env create -f workflow/envs/karyoscope.yaml
# --coverage-mode presence|graded switches to the 1 Mb-bin coverage variants.
ideogram-karyoscope: scripts/make_ideogram_karyoscope.py resources/hg19ToHg38.over.chain.gz workflow/envs/karyoscope.yaml
	conda run -n karyoscope --no-capture-output python scripts/make_ideogram_karyoscope.py

# Clean logs and temporary files
clean:
	@echo "==> Cleaning logs and temporary files..."
	rm -rf logs/
	rm -rf .snakemake/
	rm -rf results/
# 	rm -rf resources/
	rm -rf analysis/cache
	find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	@echo "==> Clean complete"

# Clean heavier local artifacts (safe to regenerate)
clean-deep: clean
	@echo "==> Cleaning deep local artifacts..."
	rm -rf .quarto/ _site/ _freeze/ _manuscript/
	rm -rf .test/
	rm -rf .pytest_cache/
	rm -f pipeline_run.html index.log
	find analysis -maxdepth 1 -type f -name "*.html" -delete 2>/dev/null || true
	find . -name ".DS_Store" -not -path "./.git/*" -delete 2>/dev/null || true
	@echo "==> Deep clean complete"
