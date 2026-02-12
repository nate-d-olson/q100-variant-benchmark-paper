# Parse-Time vs Runtime: The dip_bed Fix Explained

## The Core Problem

Snakemake evaluates code at two different times, and mixing them caused the error:

### Parse Time (DAG Construction)

- Happens when you run `snakemake -n` or `snakemake --cores 8`
- Snakemake reads all rules to understand dependencies
- Evaluates certain parameters to build the dependency graph
- **`ensure()` parameters with lambdas are evaluated here**
- Wildcard constraints haven't filtered anything yet

### Runtime (Rule Execution)

- Happens when a rule actually runs
- Only evaluates for wildcards that match the rule
- **`params` lambdas are evaluated here**
- Wildcard constraints have already filtered which wildcards match

## The Error

```python
# BEFORE: This lambda was evaluated at PARSE TIME
output:
    bed=ensure(
        "resources/benchmarksets/{benchmark}_dip.bed",
        non_empty=True,
        sha256=lambda w: config["benchmarksets"][w.benchmark]["dip_bed"]["sha256"],
        #                                                     ^^^^^^^^
        #                                                     KeyError: 'dip_bed'
        #                                                     (when benchmark=v421_grch38_smvar)
    ),
wildcard_constraints:
    benchmark="v5q_chm13_smvar|v5q_chm13_stvar|..."  # ← This can't prevent parse-time evaluation!
```

**What happened:**

1. Snakemake parses the rule
2. Evaluates `ensure()` sha256 lambda to understand the rule structure
3. Tries ALL benchmark values from config (v421, v06, AND v5q)
4. v421 and v06 don't have `dip_bed` → KeyError
5. Error occurs before wildcard constraint can filter

**Even defensive code didn't work:**

```python
sha256=lambda w: config["benchmarksets"][w.benchmark].get("dip_bed", {}).get("sha256", "")
```

The lambda itself was still being evaluated at parse time with ALL benchmark names.

## The Fix

```python
# AFTER: No lambda in ensure() - nothing evaluated at parse time
output:
    bed=ensure(
        "resources/benchmarksets/{benchmark}_dip.bed",
        non_empty=True,
        # No sha256 parameter - no parse-time evaluation!
    ),
params:
    url=lambda w: config["benchmarksets"][w.benchmark]["dip_bed"]["url"],
    sha256=lambda w: config["benchmarksets"][w.benchmark]["dip_bed"]["sha256"],
    # ↑ These lambdas evaluated at RUNTIME - after wildcard filtering!
wildcard_constraints:
    benchmark="v5q_chm13_smvar|v5q_chm13_stvar|..."  # ← Now effective!
shell:
    """
    # Manual SHA256 validation replaces ensure() validation
    echo "{params.sha256}  {output.bed}.tmp" | sha256sum -c -
    """
```

**What happens now:**

### Parse Time (safe)

1. Snakemake parses the rule
2. No lambdas in `ensure()` → nothing to evaluate
3. No error! ✅

### Runtime (safe)

1. Wildcard constraint limits matches to v5q benchmarks only
2. Rule runs with `benchmark=v5q_grch38_smvar` (for example)
3. Now `params` lambdas are evaluated
4. v5q benchmarks DO have `dip_bed` → No error! ✅

## Key Lesson

**Where you put lambdas matters in Snakemake:**

| Location | Evaluation Time | Can Use Wildcard Constraints | Safe for Missing Config Keys |
|----------|----------------|----------------------------|------------------------------|
| `ensure()` parameters | **Parse time** | ❌ No - not applied yet | ❌ No - evaluates for ALL values |
| `params` | **Runtime** | ✅ Yes - already filtered | ✅ Yes - only matching wildcards |
| `input` functions | **Runtime** | ✅ Yes - already filtered | ✅ Yes - only matching wildcards |
| `shell`/`run` | **Runtime** | ✅ Yes - already filtered | ✅ Yes - only matching wildcards |

## Visual Timeline

```
User runs: snakemake --cores 8
│
├─ PARSE TIME ───────────────────────────────────┐
│  1. Read Snakefile                             │
│  2. Parse all rules                            │
│  3. Evaluate ensure() parameters ← ERROR HERE! │
│  4. Build dependency graph (DAG)               │
│  5. Apply wildcard constraints                 │
│                                                 │
├─ RUNTIME ──────────────────────────────────────┤
│  6. Filter rules by wildcard constraints       │
│  7. Execute matching rules                     │
│  8. Evaluate params lambdas ← SAFE HERE!       │
│  9. Run shell commands                         │
└─────────────────────────────────────────────────┘
```

## Similar Snakemake Gotchas

This same issue can occur with:

- `ensure()` with any lambda parameter
- `protected()` with lambda
- `temp()` with lambda
- `ancient()` with lambda

**Solution:** Always move lambdas to `params` or `input` functions, not output wrappers.

## References

- [Snakemake documentation: Rules](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html)
- [Snakemake documentation: Wildcards](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#wildcards)
- [Snakemake documentation: Wildcard constraints](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#wildcards-constraints)
