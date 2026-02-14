# GRCh38 Debug Subset

Compact subset extracted from GRCh38 `v5.0q` benchmark resources for testing/debugging.

## Design Constraints

- One region per chromosome
- Exactly 7 chromosomes: 5 autosomes + `chrX` + `chrY`
- Includes overlap-rich regions spanning multiple exclusion/context tracks
- Includes one zero-overlap control region (`no_overlap_control`)
- Captures variant size bins from 1bp through 1-5kb with zero >5kb variants

## Regeneration

```bash
python scripts/create_grch38_debug_subset.py --force
```
