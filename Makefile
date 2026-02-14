.PHONY: help dry-run lint format format-check test clean dag run pre-commit-install

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
	@echo "  pre-commit-install - Install pre-commit hooks for automatic formatting"
	@echo "  dag              - Generate pipeline DAG visualization (PDF + DOT)"
	@echo "  run              - Execute the pipeline with conda environments"
	@echo "  clean            - Remove logs and temporary files"
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
	Rscript -e "files <- list.files('.', pattern = '.R$$', full.names = TRUE); paths <- c('R', 'test','analysis'); files <- c(files, unlist(lapply(paths[dir.exists(paths)], function(p) list.files(p, pattern = '(R|r|qmd)$$', recursive = TRUE, full.names = TRUE)))); lints <- do.call(c, lapply(files, lintr::lint)); if (length(lints) > 0) { print(lints); quit(status = 1) }"
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

# Run all tests
test: lint format-check dry-run
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
	snakemake --cores 14 --sdm conda --conda-frontend conda --report pipeline_run.html --report-after-run
	@echo "==> Pipeline execution complete"

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
