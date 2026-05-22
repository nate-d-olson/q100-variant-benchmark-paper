# Pangene Reference-Walk Layout Plan

Goal: extend `scripts/pangene_gfa_to_svg.js` with a publication-oriented layout that linearizes a Pangene GFA graph on a reference walk and highlights selected haplotypes.

## Assessment

The implementation should stay practical and avoid becoming a graph layout framework. The current GFA-to-SVG script already handles parsing and editable SVG output, so the best path is to add one new layout mode while keeping the existing Pangene-like layout unchanged.

The main uncertainty is visual readability of genes not on the reference walk. To avoid a rewrite later, isolate node placement from rendering and keep positions in a single object.

## Defaults

Keep the current default:

```text
--layout pangene
```

When the user requests:

```text
--layout reference-walk
```

default to:

```text
--reference-walk GRCh38
--highlight-walk HG002#1
--highlight-walk HG002#2
```

The common command should be:

```bash
scripts/pangene_gfa_to_svg.js \
  --input work/pangene/SULT1A1.gfa \
  --output manuscript/figs/vector/SULT1A1.hg002-linear.svg \
  --layout reference-walk \
  --title SULT1A1
```

Allow overrides for other projects:

```text
--reference-walk STR
--highlight-walk STR
```

`--highlight-walk` should be repeatable.

## Initial CLI Scope

Add:

```text
--layout pangene|reference-walk
--reference-walk STR
--highlight-walk STR
--fade-opacity FLOAT
```

Retain existing options:

```text
--block-height INT
--font-size INT
--background STR
--title STR
```

Defer:

```text
--reference-opacity
--highlight-opacity
--no-legend
--hide-unhighlighted
--layout-overrides
```

## Internal Structure

Keep placement, membership, and rendering separable:

```js
function layoutReferenceWalk(g, conf, options) {
  const walks = resolveReferenceWalks(g, options);
  const membership = computeWalkMembership(walks);
  const refLayout = layoutReferenceGenes(g, walks.reference, conf);
  const altLayout = placeOffReferenceGenes(g, walks, refLayout, conf);
  return buildRenderableGraph(refLayout, altLayout, membership);
}
```

This separation allows later manual nudges or a better alternate-gene placement heuristic without rewriting SVG output.

## Placement Heuristic

Initial deterministic heuristic:

- GRCh38/reference genes sit on a single horizontal lane.
- HG002#1-only non-reference genes go above the reference.
- HG002#2-only non-reference genes go below the reference.
- Non-reference genes present in both highlighted haplotypes go on a shared alternate lane.
- Genes absent from both highlighted haplotypes are shown faded if they are placed, but can be omitted from initial alternate placement if they are not needed to draw highlighted paths.

For each non-reference gene in a highlighted walk:

- Find the nearest previous and next genes in that walk that are present on the reference walk.
- Place the gene between those anchors.
- If multiple non-reference genes occur between the same anchors, preserve their haplotype order.
- If only one anchor exists, place the gene just before or after that anchor.

This keeps highlighted paths near the reference location where the haplotype diverges.

## Rendering

Draw:

- Faded background graph edges first.
- Highlighted haplotype paths as colored polylines through gene centers.
- Gene blocks and labels last.

Styling:

```text
HG002#1 path: #0072B2
HG002#2 path: #D55E00
genes in either highlighted haplotype: opacity 1.0
genes not in highlighted haplotypes: fade-opacity, default 0.2
```

Use stable SVG groups for later editing:

```xml
<g class="gene" data-gene="SULT1A1" data-lane="reference" opacity="1">
```

## Future Layout Feedback

If automatic placement is hard to interpret, add a manual override file later:

```bash
--layout-overrides config/pangene-layout-overrides.tsv
```

Potential format:

```text
gene    x       y       lane
SULT1A1 280     150     hg002_1
NPIPB8  228     140     hg002_shared
```

Prepare for this by keeping all positions in one `positions[sid]` object and not scattering x/y calculations across rendering code.
