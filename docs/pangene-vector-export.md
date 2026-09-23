# Pangene Gene Graph Figures (Fig 4)

Fig 4 shows two Pangene gene-graph schematics:

- **Fig 4A:** the SULT1A1 HG002 haplotype events
- **Fig 4B:** the PMS2 segmental-duplication inversion

Each panel is made in three steps:

1. Extract the local GFA subgraph from the Pangene web viewer.
2. Render it to editable SVG with `scripts/pangene_gfa_to_svg.js`.
3. Finish the layout by hand in Affinity Designer.

| Panel | Input GFA | Scripted output | Hand-finished source |
| --- | --- | --- | --- |
| Fig 4A | `data/pangene/SULT1A1.gfa` | `figures/vector/SULT1A1.events.white.{svg,pdf}` | `figures/manual/fig4a_SULT1A1.events.af` |
| Fig 4B | `data/pangene/PMS2.gfa` | `figures/vector/PMS2.events.white.{svg,pdf}` | `figures/manual/fig4b_PMS2.events.af` |

The GFAs are committed because the public Pangene server is not guaranteed to
stay available. Running the commands below on the committed GFAs reproduces
the SVGs in `figures/vector/` byte for byte (checked 2026-09-23).

## Background

Pangene's viewer (`gfa-server` from [gfatools](https://github.com/lh3/gfatools))
draws graphs on an HTML canvas, so it only produces raster images. It does
embed the local GFA subgraph in the page, in a `<textarea id="gfa-text">`
element. The two scripts below reuse that GFA:

- `scripts/extract_pangene_gfa.py` pulls the embedded GFA out of a saved
  viewer page.
- `scripts/pangene_gfa_to_svg.js` is dependency-free Node.js. It adapts the
  parser and graph-panel layout from `gfatools/js/gfa.js` and
  `gfatools/js/gfa-plot.js`, and writes SVG polygons, lines and text instead
  of drawing to a canvas.

## Reproduce Fig 4

### Step 1 (optional): re-extract the GFAs

Both GFAs come from graph `human472-1.1a2`, using the URLs cited in the Fig 4
legend: SULT1A1 with `step=3` and PMS2 with `step=10`.

```bash
for spec in SULT1A1:3 PMS2:10; do
  GENE=${spec%%:*}
  STEP=${spec##*:}
  curl -L \
    "https://pangene.bioinweb.org/view?graph=human472-1.1a2&gene=${GENE}&step=${STEP}&ori=" \
    -o "/tmp/${GENE}.pangene.html"
  scripts/extract_pangene_gfa.py \
    --input "/tmp/${GENE}.pangene.html" \
    --output "data/pangene/${GENE}.gfa"
done
```

### Step 2: render the SVGs and convert them to PDF

```bash
scripts/pangene_gfa_to_svg.js \
  --input data/pangene/SULT1A1.gfa \
  --output figures/vector/SULT1A1.events.white.svg \
  --layout event-tracks \
  --background white \
  --title "SULT1A1 HG002 haplotype events"

scripts/pangene_gfa_to_svg.js \
  --input data/pangene/PMS2.gfa \
  --output figures/vector/PMS2.events.white.svg \
  --layout pms2-events \
  --background white \
  --title "PMS2 segmental duplication inversion"

for f in SULT1A1 PMS2; do
  rsvg-convert -f pdf \
    "figures/vector/${f}.events.white.svg" \
    -o "figures/vector/${f}.events.white.pdf"
done
```

Check that `pdfinfo` reports a page width of 504 pt (7 inches).

### Step 3: finish by hand

The submitted panels were edited in Affinity Designer, starting from the
scripted SVGs. Those edits are not scripted. The `.af` files in
`figures/manual/` are the source of record for the final artwork.

## Layouts

Run `scripts/pangene_gfa_to_svg.js --help` to see all options.

- `--layout pangene`: closest to the Pangene web graph.
- `--layout reference-walk`: linearizes the graph on a reference walk
  (GRCh38 by default) and highlights the selected haplotypes
  (`--highlight-walk`, HG002#1 and HG002#2 by default).
- `--layout event-tracks`: a schematic written specifically for SULT1A1 in the
  manuscript.
  - The script checks that the GRCh38, HG002#1 and HG002#2 walks are present,
    but the event labels are hard-coded.
  - It shows the paternal SULT1A1 duplication (13 kbp insertion), the maternal
    NPIPB6+EIF3C deletion (114 kb), the two additional
    BOLA2B-SLX1A-SULT1A3 cassette copies (~102 kbp insertion in the v5.0q SV
    benchmark), and the CHM13 inversion context.
- `--layout pms2-events`: a schematic written specifically for PMS2.
  - It counts inverted walks directly from the GFA and draws the GRCh38
    17-gene block against the inverted order.
  - It annotates the PMS2 and pseudogene segmental-duplication copies.
  - Use `--highlight-walk SAMPLE#N` to pick which inverted haplotype is shown.

The event layouts follow manuscript sizing rules:

- width of 7 inches or less
- text of 10 pt or larger
- no per-gene labels inside the arrows

For a new gene, start with `pangene` or `reference-walk`. Add a new event
layout only if the figure needs collapsed flanks or custom annotations.
