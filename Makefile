.PHONY: help dry-run lint format format-check test clean dag run

QMD_FILES := $(shell find . -name '*.qmd' -not -path './.snakemake/*' -not -path './results/*' -not -path './logs/*')
MD_FILES := "**/*.md"

# Default target
help:
	@echo "Q100 Variant Benchmark Analysis - Makefile"
	@echo ""
	@echo "Available targets:"
	@echo "  help      - Show this help message"
	@echo "  dry-run   - Validate workflow with dry-run"
	@echo "  lint      - Run linting (Snakemake + R/Quarto + Markdown)"
	@echo "  format    - Format Snakemake, Quarto, and R files"
	@echo "  format-check - Check formatting without modifying files"
	@echo "  test      - Run all tests (lint + format-check + dry-run)"
	@echo "  dag       - Generate pipeline DAG visualization (PDF + DOT)"
	@echo "  run       - Execute the pipeline with conda environments"
	@echo "  clean     - Remove logs and temporary files"
	@echo ""
	@echo "Usage:"
	@echo "  make dry-run"
	@echo "  make lint"
	@echo "  make test"
	@echo "  make dag"
	@echo "  make run"

# Dry-run workflow validation
dry-run:
	@echo "==> Running workflow dry-run validation..."
	snakemake -n --quiet

# Lint workflow
lint-smk:
	@echo "==> Linting Snakemake workflow..."
	snakemake --lint
lint-r:
	@echo "==> Linting R scripts and Quarto files..."
	Rscript -e "files <- list.files('.', pattern = '.R$$', full.names = TRUE); paths <- c('R', 'test','analysis'); files <- c(files, unlist(lapply(paths[dir.exists(paths)], function(p) list.files(p, pattern = '(R|r|qmd)$$', recursive = TRUE, full.names = TRUE)))); lints <- do.call(c, lapply(files, lintr::lint)); if (length(lints) > 0) { print(lints); quit(status = 1) }"
lint-md:
	@echo "==> Linting Markdown files..."
	markdownlint $(MD_FILES)

lint: lint-smk lint-r

# Format Snakemake files
format:
	@echo "==> Formatting Snakemake files..."
	snakefmt workflow/
	@echo "==> Checking R formatting..."
	Rscript -e 'paths <- c("R","scripts","analysis"); paths <- paths[dir.exists(paths)]; if (length(paths) > 0) { lapply(paths, function(p) styler::style_dir(p, recursive = TRUE)) }'

format-check:
	@echo "==> Checking Snakefile formatting..."
	snakefmt --check workflow/
	@echo "==> Checking R formatting..."
	Rscript -e 'paths <- c("R","scripts"); paths <- paths[dir.exists(paths)]; if (length(paths) > 0) { lapply(paths, function(p) styler::style_dir(p, recursive = TRUE, dry = "fail")) }'

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
	snakemake --cores 14 --sdm conda --conda-frontend conda
	@echo "==> Pipeline execution complete"

# Clean logs and temporary files
clean:
	@echo "==> Cleaning logs and temporary files..."
	rm -rf logs/
	rm -rf .snakemake/
	rm -rf results/
	rm -rf resources/
	rm -rf analysis/cache
	find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	@echo "==> Clean complete"
