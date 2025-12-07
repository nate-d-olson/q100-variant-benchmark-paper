.PHONY: help dry-run lint format test clean dag

# Default target
help:
	@echo "Q100 Variant Benchmark Analysis - Makefile"
	@echo ""
	@echo "Available targets:"
	@echo "  help      - Show this help message"
	@echo "  dry-run   - Validate workflow with dry-run"
	@echo "  lint      - Run workflow linting"
	@echo "  format    - Format Snakemake files with snakefmt"
	@echo "  test      - Run all tests (lint + format check + dry-run)"
	@echo "  dag       - Generate pipeline DAG visualization (PDF + DOT)"
	@echo "  clean     - Remove logs and temporary files"
	@echo ""
	@echo "Usage:"
	@echo "  make dry-run"
	@echo "  make lint"
	@echo "  make test"
	@echo "  make dag"

# Dry-run workflow validation
dry-run:
	@echo "==> Running workflow dry-run validation..."
	snakemake -n --quiet

# Lint workflow
lint:
	@echo "==> Linting Snakemake workflow..."
	snakemake --lint
	@echo "==> Checking Snakefile formatting..."
	snakefmt --check workflow/

# Format Snakemake files
format:
	@echo "==> Formatting Snakemake files..."
	snakefmt workflow/

# Run all tests
test: lint dry-run
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

# Clean logs and temporary files
clean:
	@echo "==> Cleaning logs and temporary files..."
	rm -rf logs/
	rm -rf .snakemake/
	find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	@echo "==> Clean complete"
