# Data Loading Improvement Proposals

## Overview

This directory contains comprehensive proposals for improving data loading in Quarto analysis notebooks for the q100-variant-benchmark-paper project. The proposals were developed to address the need to make it easier to finalize figures and tables for the scientific manuscript describing the new v5.0q benchmark sets.

## Contents

### Proposals

1. **[proposal_1_caching_library.md](./proposal_1_caching_library.md)** (9 KB)
   - **Title:** Enhanced Data Loading Library with Intelligent Caching
   - **Key Idea:** Add transparent disk-based caching to existing data loading functions
   - **Performance:** 80-95% faster re-runs
   - **Implementation:** 1-2 weeks
   - **Breaking Changes:** None
   - **Status:** ✅ **SELECTED**

2. **[proposal_2_data_access_layer.md](./proposal_2_data_access_layer.md)** (15 KB)
   - **Title:** Standardized Data Access Layer with Validation Framework
   - **Key Idea:** Create comprehensive validation with explicit schemas
   - **Focus:** Data quality assurance, type safety
   - **Implementation:** 3-4 weeks
   - **Breaking Changes:** Optional migration
   - **Status:** Future consideration

3. **[proposal_3_data_snapshots.md](./proposal_3_data_snapshots.md)** (18 KB)
   - **Title:** Analysis-Ready Data Snapshots with Version Control
   - **Key Idea:** Generate versioned, pre-processed data artifacts
   - **Performance:** <1 second loads (instant)
   - **Implementation:** 2-3 weeks
   - **Breaking Changes:** Requires workflow change
   - **Status:** Future consideration (Phase 2)

### Analysis & Selection

1. **[proposal_comparison_and_selection.md](./proposal_comparison_and_selection.md)** (12 KB)
   - Detailed comparison across 6 criteria (implementation, performance, UX, maintainability, reproducibility, alignment)
   - Weighted scoring methodology (Proposal 1: 91.25, Proposal 2: 71.25, Proposal 3: 78.5)
   - Rationale for selecting Proposal 1
   - Future roadmap for other proposals

### Implementation

1. **[implementation_plan_caching.md](./implementation_plan_caching.md)** (26 KB)
   - Complete 4-phase implementation plan for Proposal 1
   - Timeline: 10-14 days
   - Detailed tasks with acceptance criteria
   - Code samples, test plans, and documentation templates
   - Risk assessment and rollback strategy
   - Success metrics and validation procedures

2. **[cache_format_comparison.md](./cache_format_comparison.md)** (20 KB) ⭐ NEW
   - **Response to PR comment:** Comparison of cache formats (RDS vs Arrow/Parquet vs SQLite/DuckDB)
   - Performance benchmarks for different data sizes
   - Memory usage analysis for exploratory analysis
   - **Recommendation:** Arrow/Parquet for 2-4x performance, optional DuckDB for queries
   - Use case analysis and migration path

3. **[implementation_plan_parquet_addendum.md](./implementation_plan_parquet_addendum.md)** (19 KB) ⭐ NEW
   - **Extension to implementation plan:** Adds Arrow/Parquet format support
   - Multi-format caching (RDS + Parquet) with auto-selection
   - Optional DuckDB query support for large datasets
   - Maintains 100% backward compatibility with RDS-only design
   - Graceful fallback when arrow package not installed

## Selection Summary

**Selected Proposal: #1 - Enhanced Data Loading Library with Intelligent Caching**

### Why Proposal 1?

1. **Perfect alignment with current phase** - Project is in manuscript finalization with frequent notebook re-runs
2. **Zero disruption** - Existing notebooks work unchanged with transparent caching
3. **Immediate impact** - 80-95% performance improvement in 1-2 weeks
4. **Low risk** - Small code changes, easy to remove if issues arise
5. **Flexibility preserved** - Maintains full exploratory analysis capabilities
6. **Foundation for future** - Can be combined with other proposals later

### Key Benefits

- ✅ **Performance:** Cached data loads in milliseconds vs. seconds
- ✅ **Backward Compatible:** Zero changes to existing notebook code
- ✅ **Automatic Invalidation:** Cache refreshes when source files change
- ✅ **Simple Management:** `clear_analysis_cache()` and `cache_stats()` utilities
- ✅ **Transparent:** Developers barely notice caching is happening
- ✅ **Low Maintenance:** ~200 lines of code, single file changes

## Cache Format Enhancement (Updated 2025-02-07)

Following PR feedback, the implementation plan has been enhanced to support multiple cache formats:

### Format Options

| Format | Performance | Compression | Memory | Best For |
|--------|-------------|-------------|--------|----------|
| **RDS** | Fast | Good | High | Small data (<5 MB) |
| **Parquet** | 2-4x faster | Excellent | Low | Medium/large data |
| **DuckDB** | Query-based | Good | Very Low | Exploratory queries |

### Key Enhancements

1. **Multi-format support:** RDS (default) + Parquet (optional, requires arrow package)
2. **Auto-selection:** Automatically chooses format based on data size
3. **Backward compatible:** Works without arrow, falls back to RDS gracefully
4. **Optional queries:** DuckDB integration for SQL-based exploratory analysis

### Performance Benefits

**Parquet vs RDS:**

- Read speed: 2-4x faster
- Compression: 40-60% smaller files
- Memory: Column subsetting reduces memory usage
- Queries: Direct SQL queries without full data load (with DuckDB)

**See:**

- [cache_format_comparison.md](./cache_format_comparison.md) - Detailed analysis
- [implementation_plan_parquet_addendum.md](./implementation_plan_parquet_addendum.md) - Implementation details

---

## Current Status

- **Date Created:** 2025-02-07
- **Created By:** GitHub Copilot Agent
- **Branch:** `copilot/improve-data-loading-methods`
- **Commit:** `17b956e`
- **Status:** ✅ Proposals complete, ready for implementation

## Next Steps

1. **Review proposals** with project team
2. **Approve implementation plan** for Proposal 1
3. **Create feature branch** (`feat/caching-library`)
4. **Begin Phase 1** of implementation plan
5. **Track progress** through daily standups
6. **Demo after Phase 2** (function integration)

## Implementation Timeline

| Phase | Duration | Deliverables |
|-------|----------|--------------|
| Phase 1: Core Infrastructure | Days 1-3 | Caching utilities, cache management functions |
| Phase 2: Function Integration | Days 4-7 | 4 major functions with caching enabled |
| Phase 3: Documentation & Testing | Days 8-10 | Complete docs, comprehensive test suite |
| Phase 4: Validation & Rollout | Days 11-14 | Performance benchmarks, migration guide |

**Total Estimated Time:** 10-14 days

## Context

### Project Phase

The q100-variant-benchmark-paper project is currently in the manuscript finalization phase where:

- Figures are being iteratively refined based on feedback
- Notebooks are being re-run frequently (5-10+ times per day)
- External collaborators need to review results
- Pipeline outputs are relatively stable

### Current Pain Points

1. Loading data from pipeline outputs is slow (2-5 seconds per dataset)
2. No caching mechanism exists, requiring full reload every time
3. Some notebooks use inconsistent data loading patterns
4. Large variant tables (GB-scale) are very slow to load

### Goal

Make it easier to finalize the initial set of figures and tables for the scientific manuscript while maintaining flexibility for exploratory analysis.

## Technical Details

### Current Architecture

- Data loading functions in `R/data_loading.R`
- Snakemake pipeline outputs in `results/` directory
- Two-tier approach: aggregated metrics (fast) vs. full variant tables (slow)
- Helper functions provide primary interface

### Proposed Enhancement (Proposal 1)

- Add internal `.cache_wrapper()` function
- Wrap existing loaders with caching layer
- Store cached objects in `analysis/cache/` (gitignored)
- Use content-based cache keys (file mtime + size)
- Provide management utilities (`clear_analysis_cache()`, `cache_stats()`)

### Performance Expectations

- First load: 2-5 seconds (same as current)
- Cached load: 0.1-0.5 seconds (10-50x faster)
- Overall improvement: 80-95% for typical notebook re-runs

## References

### Related Documentation

- [Pipeline Outputs Reference](../pipeline-outputs.md) - Structure of Snakemake outputs
- [Data Dictionary](../data-dictionary.md) - Metric definitions
- [API Reference](../api-reference.md) - Function documentation
- [Architecture](../architecture.md) - System design overview

### Code Files

- `R/data_loading.R` - Current data loading functions (to be enhanced)
- `analysis/*.qmd` - Quarto notebooks using data loading functions
- `config/config.yaml` - Pipeline configuration

## Questions?

For questions or clarifications about these proposals:

1. Review the detailed proposal documents
2. Check the implementation plan for technical details
3. Consult the comparison document for selection rationale
4. Refer to the project README for general context

---

**Document Version:** 1.0  
**Last Updated:** 2025-02-07  
**Author:** GitHub Copilot Agent  
**Status:** Complete and ready for implementation
