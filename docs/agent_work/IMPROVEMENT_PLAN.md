# Q100 Variant Benchmark Pipeline - Improvement Plan

**Date:** 2026-01-16
**Status:** Proposed
**Priority:** High Impact, High Value

## Overview

This document outlines 5 key improvements to enhance the maintainability, reliability, and quality of the Q100 variant benchmark pipeline codebase.

---

## Feature 1: Code Duplication Removal

### Problem

`workflow/rules/common.smk` contains duplicate function definitions (lines 29-69 duplicated at lines 86-129):

- `get_comparison_files()` - appears twice
- `get_stratifications_for_comp()` - appears twice
- `get_strat_inputs()` - appears twice with minor variations

This creates maintenance burden and risk of inconsistencies when updating logic.

### Solution

- Remove duplicate function definitions
- Consolidate into single canonical implementations
- Add clear documentation for each function's purpose and parameters
- Add type hints for better IDE support

### Benefits

- **Maintainability**: Single source of truth for helper functions
- **Consistency**: No risk of divergent implementations
- **Readability**: Clearer code structure with ~50 fewer lines

### Implementation Scope

- File: `workflow/rules/common.smk`
- Estimated lines removed: ~40
- Complexity: Low
- Risk: Low (functions are identical)

---

## Feature 2: Python Code Quality Tooling

### Problem

The codebase lacks automated Python code quality enforcement:

- No linting (ruff not configured)
- No auto-formatting (black not configured)
- No type checking (mypy not used despite type hints present)
- Inconsistent code style across Python scripts

### Solution

Implement comprehensive Python tooling suite:

1. **Ruff for Linting**
   - Fast, modern Python linter
   - Replace multiple tools (flake8, isort, etc.)
   - Configure in `pyproject.toml`

2. **Black for Formatting**
   - Consistent code formatting
   - Integrate with Makefile and pre-commit

3. **mypy for Type Checking**
   - Static type verification
   - Catch type errors before runtime
   - Enhance IDE support

4. **Pre-commit Hooks**
   - Automatic checks before commits
   - Prevent bad code from entering repository

### Configuration Files to Add

```toml
# pyproject.toml additions
[tool.ruff]
line-length = 100
target-version = "py39"
select = ["E", "F", "I", "N", "W", "B", "C90"]

[tool.black]
line-length = 100
target-version = ['py39']

[tool.mypy]
python_version = "3.9"
warn_return_any = true
warn_unused_configs = true
```

### Benefits

- **Code Quality**: Catch bugs and style issues early
- **Consistency**: Uniform code style across project
- **Developer Experience**: Better IDE support with type hints
- **Collaboration**: Easier code reviews with automated checks

### Implementation Scope

- Add configuration files: `pyproject.toml`, `.pre-commit-config.yaml`
- Update `Makefile` with new targets: `make lint-py`, `make format-py`, `make typecheck`
- Fix existing linting/formatting issues
- Complexity: Medium
- Risk: Low (non-breaking changes)

---

## Feature 3: CI/CD Pipeline with GitHub Actions

### Problem

No continuous integration or automated testing:

- Tests not run automatically on commits/PRs
- No validation of pull requests before merge
- Manual testing burden on developers
- Risk of breaking changes reaching main branch

### Solution

Implement GitHub Actions CI/CD pipeline:

**Workflows to Create:**

1. **CI Workflow** (`.github/workflows/ci.yml`)
   - Trigger: On push and pull request
   - Jobs:
     - Lint Python code (ruff)
     - Format check (black --check)
     - Type check (mypy)
     - Run unit tests (pytest)
     - Validate Snakemake workflow
     - Check Snakemake formatting (snakefmt --check)

2. **Test Coverage Workflow**
   - Generate coverage reports
   - Upload to Codecov or similar service
   - Comment coverage changes on PRs

3. **Documentation Build** (optional)
   - Build and validate documentation
   - Deploy to GitHub Pages

**Example CI Job Structure:**

```yaml
name: CI

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: '3.11'
      - name: Install dependencies
        run: |
          pip install -r requirements.txt
          pip install pytest pytest-cov ruff black mypy
      - name: Lint
        run: ruff check .
      - name: Format check
        run: black --check .
      - name: Type check
        run: mypy workflow/scripts/
      - name: Run tests
        run: pytest --cov=workflow --cov-report=xml
```

### Benefits

- **Quality Assurance**: Automated testing catches issues early
- **Confidence**: Safe to merge PRs with passing checks
- **Documentation**: CI badge shows build status
- **Efficiency**: Reduces manual testing burden

### Implementation Scope

- Create `.github/workflows/ci.yml`
- Update README with CI badge
- Configure test coverage reporting
- Complexity: Medium
- Risk: Low (doesn't affect existing code)

---

## Feature 4: Comprehensive Test Coverage

### Problem

Testing coverage is minimal:

- Only validators and exceptions have tests
- **0%** coverage of Python scripts (14 scripts, 1,836 lines)
- **0%** integration tests for Snakemake rules
- **0%** tests for data processing logic
- Limited test fixtures (only 2: sample.bed, sample.vcf)

### Solution

Expand testing in phases:

**Phase 1: Unit Tests for Python Scripts**

- `count_variants_by_type.py` - test variant counting logic
- `stratify_comparison.py` - test stratification processing
- `combine_beds_with_id.py` - test BED interval merging
- `expand_annotations.py` - test annotation expansion
- Target: 80% line coverage for scripts/

**Phase 2: Integration Tests**

- Test key Snakemake rules with realistic data
- Rules to test:
  - `download_vcf` - test download and checksum validation
  - `materialize_exclusion` - test BED file processing
  - `compute_exclusion_metrics` - test metric calculation
  - `aggregate_comparison` - test result aggregation

**Phase 3: Expanded Test Fixtures**

- Add multi-reference test data (GRCh37, GRCh38, CHM13v2.0)
- Add edge case fixtures (empty files, malformed data)
- Add realistic benchmark set samples

**Test Structure:**

```
tests/
├── unit/
│   ├── test_validators.py (existing)
│   ├── test_count_variants.py (new)
│   ├── test_stratify.py (new)
│   └── test_combine_beds.py (new)
├── integration/
│   ├── test_download_rules.py (new)
│   ├── test_exclusion_rules.py (new)
│   └── test_comparison_rules.py (new)
└── fixtures/
    ├── vcf/ (new)
    ├── bed/ (new)
    └── config/ (new)
```

### Benefits

- **Reliability**: Catch bugs before production
- **Refactoring Safety**: Confidently modify code
- **Documentation**: Tests show expected behavior
- **Quality Metrics**: Track coverage trends

### Implementation Scope

- Add 15-20 new test files
- Create comprehensive test fixtures
- Update CI to run all tests
- Complexity: High
- Risk: Low (tests don't affect production code)

---

## Feature 5: Architecture & API Documentation

### Problem

Documentation gaps prevent easy understanding and contribution:

- No architecture overview or data flow diagrams
- No API reference for Python scripts
- Configuration schema undocumented
- No troubleshooting guide for common failures
- No development setup instructions

### Solution

Create comprehensive documentation suite:

**1. Architecture Documentation** (`docs/architecture.md`)

- High-level pipeline overview diagram
- Data flow through rules
- Directory structure explanation
- Rule dependency graph
- Resource requirements (CPU, memory, disk)

**2. API Reference** (`docs/api-reference.md`)

- Automatically generated from docstrings
- Document all Python scripts and functions
- Include type signatures and examples
- Cross-reference with Snakemake rules

**3. Configuration Guide** (`docs/configuration.md`)

- Explain config.yaml structure
- Document all config parameters
- Provide examples for common use cases:
  - Adding new benchmark sets
  - Adding new stratifications
  - Customizing comparisons
- Document config.schema.yaml validation rules

**4. Troubleshooting Guide** (`docs/troubleshooting.md`)

- Common errors and solutions
- Debugging workflow failures
- Performance optimization tips
- FAQ section

**5. Developer Guide** (`docs/development.md`)

- Development environment setup
- Code style guidelines
- Testing best practices
- Contributing guidelines
- Release process

**6. Enhanced README**

- Quick start guide
- Installation instructions
- Basic usage examples
- Link to detailed documentation

### Documentation Tools

- Use **mkdocs** or **Sphinx** for documentation site
- Auto-generate API docs from docstrings
- Include Mermaid diagrams for workflows
- Host on GitHub Pages

### Benefits

- **Onboarding**: New contributors get up to speed faster
- **Usability**: Users understand how to customize pipeline
- **Maintenance**: Clear architecture aids debugging
- **Collaboration**: Better knowledge sharing

### Implementation Scope

- Create 6 new documentation files
- Set up documentation build system
- Generate API reference
- Create workflow diagrams
- Complexity: High
- Risk: Low (documentation only)

---

## Implementation Roadmap

### Phase 1: Quick Wins (Week 1)

1. **Feature 1**: Remove code duplication (~2 hours)
2. **Feature 2**: Add Python tooling configuration (~4 hours)

### Phase 2: Automation (Week 2)

1. **Feature 3**: Implement CI/CD pipeline (~8 hours)

### Phase 3: Testing (Weeks 3-4)

1. **Feature 4**: Expand test coverage (~16 hours)
   - Week 3: Unit tests for scripts
   - Week 4: Integration tests

### Phase 4: Documentation (Week 5)

1. **Feature 5**: Create comprehensive documentation (~12 hours)

**Total Estimated Effort:** ~42 hours over 5 weeks

---

## Success Metrics

### Code Quality

- **Duplication**: 0 duplicate functions in common.smk
- **Linting**: 100% of Python files pass ruff checks
- **Formatting**: 100% of Python files formatted with black
- **Type Coverage**: 100% of Python files type-checked with mypy

### Testing

- **Unit Test Coverage**: ≥80% for workflow/scripts/
- **Integration Tests**: ≥5 critical rules tested
- **CI Pass Rate**: 100% of main branch commits pass CI

### Documentation

- **Coverage**: All 6 documentation sections complete
- **API Docs**: 100% of public functions documented
- **Examples**: ≥3 configuration examples provided

---

## Risks & Mitigation

| Risk | Impact | Probability | Mitigation |
|------|--------|-------------|------------|
| Breaking changes during refactoring | High | Low | Comprehensive testing before merge |
| CI/CD slows down development | Medium | Low | Optimize CI runtime, use caching |
| Documentation becomes outdated | Medium | Medium | Automate API doc generation |
| Testing effort exceeds estimate | Low | Medium | Prioritize critical paths first |

---

## Next Steps

1. **Review**: Team review of this improvement plan
2. **Prioritization**: Confirm priority order of features
3. **Branch Creation**: Create `feature/codebase-improvements` branch
4. **Implementation**: Begin with Phase 1 (Quick Wins)
5. **Iteration**: Regular check-ins and adjustments

---

## References

- Original codebase analysis: See exploration agent output
- Snakemake best practices: <https://snakemake.readthedocs.io/>
- Python tooling: ruff, black, mypy documentation
- Testing: pytest best practices
