Improve and simplify vcf processing using functions from the truvari api to replace code in python scripts when appropriate.

# Project To Dos

- Update variant cache schema

## Analysis Notebooks

- update input and outputs for notebooks are consistent with code

### benchmark_difficult.qmd

- Revise diff coverage plots and fix inputs (coverage data loading/plotting pipeline)
- Check X7 as table column — unexpected column in diff_cov_df after loading
- Check warning message: "expected pieces with additional discarded pieces for 24 rows" in separate() call — confirm this is all rows and fix parsing
- Revise dead code fold-change chunks (eval=FALSE) using `strat_ids` — starting point for new figures

### benchmark_exclusions.qmd

- Add upset plot for exclusion region intersections

### benchmarkset_characterization.qmd

- Add variants to the exclusion tables
- Small variant counts breakdown by <15bp , 15 - 49bp
- fix PP missing from benchmark region size distribution fig
- fix difficult regions included and difficult region bases excluded code

Documentation
- Column descriptions for variant tables

Analyses
- Look for small intervals in v5q

Low priority
- Optimize variant tables: data.table, database, arrow, etc.?
