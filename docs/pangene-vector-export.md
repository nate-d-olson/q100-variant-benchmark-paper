# Pangene Gene Graph Vector Export Options

Goal: recreate the Pangene web viewer graph for `SULT1A1` as an editable SVG or PDF for publication figures.

Target web view:

```text
https://pangene.bioinweb.org/view?graph=human472-1.1a2&gene=SULT1A1&step=3&ori=
```

Current local raster reference:

```text
manuscript/figs/gene1.png
```

## What I Verified

The Pangene source repo documents graph visualization through `gfa-server`, which is distributed with `gfatools` and bundled in the Pangene binary release. The relevant code paths are:

- `pangene/README.md`: visualization uses `gfa-server -d html data/*.gfa.gz`.
- `gfatools/gfa-server.go`: `/view` extracts a local GFA subgraph with `gfa_extract()`, embeds the result in a readonly `<textarea id="gfa-text">`, and calls JavaScript plotting functions.
- `gfatools/js/gfa.js`: parses the embedded GFA into `seg`, `arc`, and `walk` objects.
- `gfatools/js/gfa-plot.js`: computes graph layout and draws with `canvas.getContext("2d")`.

The viewer is raster-only because it draws into HTML canvas. However, the plot itself is made from vector-like primitives: line segments, arrow polygons, text labels, colors, and fixed coordinates. That makes native SVG export feasible without changing Pangene graph construction.

## Recommendation

For the manuscript, the best route is Approach 2: make a small standalone exporter that reuses the existing Pangene/GFA layout logic and writes SVG. It keeps the webapp untouched, gives a reproducible command for the exact gene figure, and avoids browser-specific canvas export.

Approach 1 is best if you want to contribute a general feature upstream to `gfatools`. Approach 3 is the fastest fallback, but it will not exactly match Pangene's layout unless tuned.

## Approach 1: Add Native SVG Export To The Web Viewer

This modifies `gfatools/js/gfa-plot.js` so the same page can download SVG. The existing renderer already computes all positions in `gfa_plot_cal_pos()`, colors in `gfa_plot_rank2color()`, and arrow dimensions in `gfa_plot_arrow()`. Add parallel SVG functions instead of trying to convert the canvas bitmap.

Suggested workflow:

```bash
git clone https://github.com/lh3/gfatools.git
cd gfatools

rg -n "gfa_plot_graph|gfa_plot_arrow|canvas|getContext|gfa_plot_cal_pos" js/gfa-plot.js gfa-server.go
```

Implementation shape:

```js
function gfa_svg_arrow(svg, x, y, len, w, rev, text, fs, color_stroke, color_fill) {
  const pts = rev
    ? [[x, y], [x + w, y - w], [x + w + len, y - w], [x + len, y], [x + w + len, y + w], [x + w, y + w]]
    : [[x, y], [x - w, y - w], [x - w + len, y - w], [x + len, y], [x - w + len, y + w], [x - w, y + w]];

  const poly = document.createElementNS("http://www.w3.org/2000/svg", "polygon");
  poly.setAttribute("points", pts.map(p => p.join(",")).join(" "));
  if (color_fill) poly.setAttribute("fill", color_fill);
  else poly.setAttribute("fill", "none");
  if (color_stroke) poly.setAttribute("stroke", color_stroke);
  svg.appendChild(poly);

  if (text != null) {
    const t = document.createElementNS("http://www.w3.org/2000/svg", "text");
    t.setAttribute("x", x + len / 2);
    t.setAttribute("y", y - w - 2);
    t.setAttribute("text-anchor", "middle");
    t.setAttribute("font-family", "monospace");
    t.setAttribute("font-size", fs || 9);
    t.setAttribute("fill", color_stroke || color_fill || "#000000");
    t.textContent = text;
    svg.appendChild(t);
  }
}
```

Then add a `gfa_plot_graph_svg(conf, g)` function that mirrors `gfa_plot_graph()`:

```js
const ret = gfa_plot_cal_pos(conf, g);
const sub = ret[0], pos = ret[1];
const r2c = gfa_plot_rank2color(g);
```

Draw SVG `<line>` elements for edges and SVG `<polygon>` plus `<text>` for genes. Add a download button in `gfa-server.go` near the existing `Replot` button.

Example local test commands:

```bash
make kalloc.o gfa-base.o gfa-io.o gfa-util.o
go build -o gfa-server gfa-server.go

./gfa-server \
  -j js \
  -p 8000 \
  /path/to/human472-1.1a2.gfa.gz

open "http://127.0.0.1:8000/view?graph=human472-1.1a2&gene=SULT1A1&step=3&ori="
```

Export PDF from SVG:

```bash
inkscape SULT1A1.pangene.svg --export-filename=SULT1A1.pangene.pdf
```

Pros:

- Closest to the existing webapp.
- Could support both graph and walk panels.
- Useful upstream feature.

Cons:

- Requires editing JavaScript and Go template output.
- Still browser-centered unless you also add a command-line path.

## Approach 2: Standalone Command-Line SVG Exporter

This is the simplest reproducible manuscript workflow. Use `gfa-server` or the public Pangene server to get the exact local GFA subgraph, then render it to SVG with a script that ports the small amount of drawing logic from `gfa-plot.js`.

The implemented workflow has two small scripts:

- `scripts/extract_pangene_gfa.py`: extracts the local GFA subgraph embedded in the Pangene viewer HTML.
- `scripts/pangene_gfa_to_svg.js`: renders the graph panel to editable SVG using Pangene/GFA viewer layout logic.

Create an output directory:

```bash
mkdir -p work/pangene manuscript/figs/vector
```

Fetch the rendered page from the public server:

```bash
curl -L \
  "https://pangene.bioinweb.org/view?graph=human472-1.1a2&gene=SULT1A1&step=3&ori=" \
  -o work/pangene/SULT1A1.pangene.html
```

Extract the embedded GFA from the textarea:

```bash
scripts/extract_pangene_gfa.py \
  --input work/pangene/SULT1A1.pangene.html \
  --output work/pangene/SULT1A1.gfa
```

Render SVG:

```bash
scripts/pangene_gfa_to_svg.js \
  --input work/pangene/SULT1A1.gfa \
  --output manuscript/figs/vector/SULT1A1.pangene.svg \
  --title SULT1A1
```

The renderer in `scripts/pangene_gfa_to_svg.js` is dependency-free Node.js. It adapts only the parser and graph-panel layout logic needed from `gfatools/js/gfa.js` and `gfatools/js/gfa-plot.js`, then emits SVG `<line>`, `<polygon>`, and `<text>` elements instead of drawing to canvas.

Gene blocks are written as explicit SVG polygons with hex colors, matching strokes, and non-scaling stroke widths. This is more robust in vector editors than relying on browser-only canvas output or CSS `hsl(...)` colors.

Useful options:

```bash
scripts/pangene_gfa_to_svg.js --help
```

Common examples:

```bash
# Transparent background, best as an editable source figure.
scripts/pangene_gfa_to_svg.js \
  --input work/pangene/SULT1A1.gfa \
  --output manuscript/figs/vector/SULT1A1.pangene.svg \
  --title SULT1A1

# White background, useful for quick previewing.
scripts/pangene_gfa_to_svg.js \
  --input work/pangene/SULT1A1.gfa \
  --output manuscript/figs/vector/SULT1A1.pangene.white.svg \
  --background white \
  --title SULT1A1

# Slightly taller gene blocks if the default looks too thin in a journal figure.
scripts/pangene_gfa_to_svg.js \
  --input work/pangene/SULT1A1.gfa \
  --output manuscript/figs/vector/SULT1A1.pangene.tall.svg \
  --block-height 14 \
  --title SULT1A1
```

Reference-walk layout:

```bash
# Linearize on GRCh38 and highlight HG002#1/HG002#2.
# GRCh38, HG002#1, and HG002#2 are defaults when --layout reference-walk is used.
scripts/pangene_gfa_to_svg.js \
  --input work/pangene/SULT1A1.gfa \
  --output manuscript/figs/vector/SULT1A1.hg002-linear.svg \
  --layout reference-walk \
  --title SULT1A1

# White background variant for previewing.
scripts/pangene_gfa_to_svg.js \
  --input work/pangene/SULT1A1.gfa \
  --output manuscript/figs/vector/SULT1A1.hg002-linear.white.svg \
  --layout reference-walk \
  --background white \
  --title SULT1A1
```

Reference-walk mode arranges reference genes in GRCh38 walk order, draws HG002#1 and HG002#2 as highlighted paths, and fades genes that are absent from the highlighted walks using `--fade-opacity` when such genes are present. You can override the defaults:

```bash
scripts/pangene_gfa_to_svg.js \
  --input work/pangene/MYGENE.gfa \
  --output manuscript/figs/vector/MYGENE.linear.svg \
  --layout reference-walk \
  --reference-walk CHM13#0 \
  --highlight-walk SAMPLE#1 \
  --highlight-walk SAMPLE#2 \
  --fade-opacity 0.15 \
  --title MYGENE
```

The off-reference placement heuristic is deliberately simple: genes not matched to the reference walk are placed near their nearest reference-walk anchors in the highlighted haplotype path, above or below the reference lane. If a future figure is hard to interpret, the next planned extension is a small manual layout override table rather than replacing the parser or SVG renderer.

Manuscript event-track layout:

The event-track layouts are formatted for manuscript placement: SVG physical width is capped at 7 inches, visible text uses 10 pt or larger font sizes, individual gene names inside arrows are omitted, and gene blocks are drawn taller than the Pangene-like graph export.

```bash
# Compact schematic for the SULT1A1 HG002 haplotype events described in the manuscript.
scripts/pangene_gfa_to_svg.js \
  --input work/pangene/SULT1A1.gfa \
  --output manuscript/figs/vector/SULT1A1.events.svg \
  --layout event-tracks \
  --title "SULT1A1 HG002 haplotype events"

# White-background variant, useful for previewing and direct manuscript placement.
scripts/pangene_gfa_to_svg.js \
  --input work/pangene/SULT1A1.gfa \
  --output manuscript/figs/vector/SULT1A1.events.white.svg \
  --layout event-tracks \
  --background white \
  --title "SULT1A1 HG002 haplotype events"
```

Convert the event schematic to PDF:

```bash
rsvg-convert -f pdf \
  manuscript/figs/vector/SULT1A1.events.svg \
  -o manuscript/figs/vector/SULT1A1.events.pdf

rsvg-convert -f pdf \
  manuscript/figs/vector/SULT1A1.events.white.svg \
  -o manuscript/figs/vector/SULT1A1.events.white.pdf
```

The `event-tracks` layout is intentionally manuscript-specific. It verifies that the input GFA contains the expected GRCh38, HG002#1, and HG002#2 walks and uses gene colors from the GFA, but it draws a compact schematic of the described events rather than trying to infer every event label generically. The output emphasizes:

- GRCh38, HG002#1, and HG002#2 as separate tracks.
- Collapsed syntenic blocks for flanking regions that are not central to the figure.
- The HG002#1 paternal `SULT1A1` duplication labeled as a 13 kbp insertion.
- The HG002#2 maternal `NPIPB6` plus `EIF3C` deletion labeled as 114 kb.
- Two additional `BOLA2B-SLX1A-SULT1A3` cassette copies in each HG002 haplotype relative to GRCh38, labeled as the approximately 102 kbp insertion in the v5.0q SV benchmark.
- The shared CHM13 inversion annotation explaining why that region is excluded from the CHM13-based benchmark.

PMS2 inversion event layout:

```bash
# Compact schematic for the PMS2 segmental-duplication-mediated inversion.
scripts/pangene_gfa_to_svg.js \
  --input pms2.gfa \
  --output manuscript/figs/vector/PMS2.events.svg \
  --layout pms2-events \
  --title "PMS2 segmental duplication inversion"

# White-background variant for previewing and direct manuscript placement.
scripts/pangene_gfa_to_svg.js \
  --input pms2.gfa \
  --output manuscript/figs/vector/PMS2.events.white.svg \
  --layout pms2-events \
  --background white \
  --title "PMS2 segmental duplication inversion"
```

Convert the PMS2 schematic to PDF:

```bash
rsvg-convert -f pdf \
  manuscript/figs/vector/PMS2.events.svg \
  -o manuscript/figs/vector/PMS2.events.pdf

rsvg-convert -f pdf \
  manuscript/figs/vector/PMS2.events.white.svg \
  -o manuscript/figs/vector/PMS2.events.white.pdf
```

The `pms2-events` layout is a PMS2-specific schematic. It counts the inverted walks directly from `pms2.gfa`, chooses a representative inverted haplotype, draws the GRCh38 17-gene block against the inverted order, and annotates the duplicated PMS2/pseudogene segmental-duplication copies that can cause minimap2 to swap mappings between PMS2 and its pseudogene. To force a specific inverted walk instead of the first detected one, pass it explicitly:

```bash
scripts/pangene_gfa_to_svg.js \
  --input pms2.gfa \
  --output manuscript/figs/vector/PMS2.events.HG00735.svg \
  --layout pms2-events \
  --highlight-walk HG00735#1 \
  --title "PMS2 segmental duplication inversion"
```

### Choosing A Layout

Use the simplest layout that answers the figure question:

- `--layout pangene`: closest to the Pangene web graph. Best for an editable vector version of the webapp view.
- `--layout reference-walk`: linearizes on a reference walk, defaults to GRCh38, and highlights selected sample haplotypes. Best when the biological point is about how one or two haplotypes differ from a reference order.
- `--layout event-tracks`: SULT1A1-specific manuscript schematic. Best for the current SULT1A1 figure because it encodes the described HG002 insertion, duplication, deletion, and CHM13 inversion context.
- `--layout pms2-events`: PMS2-specific manuscript schematic. Best for the PMS2 segmental-duplication inversion figure because it counts inverted walks from the GFA and shows the 17-gene inversion compactly.

For new genes, start with `--layout pangene` and `--layout reference-walk`. Only add a new event layout if the figure needs collapsed syntenic regions, custom event labels, or manuscript-specific interpretation that should not be inferred automatically from arbitrary GFAs.

### Reusable Commands

General Pangene-style vector export from a local GFA:

```bash
GENE=MYGENE

scripts/pangene_gfa_to_svg.js \
  --input "work/pangene/${GENE}.gfa" \
  --output "manuscript/figs/vector/${GENE}.pangene.svg" \
  --title "${GENE}"

rsvg-convert -f pdf \
  "manuscript/figs/vector/${GENE}.pangene.svg" \
  -o "manuscript/figs/vector/${GENE}.pangene.pdf"
```

Reference-walk figure with custom haplotypes:

```bash
GENE=MYGENE

scripts/pangene_gfa_to_svg.js \
  --input "work/pangene/${GENE}.gfa" \
  --output "manuscript/figs/vector/${GENE}.reference-walk.svg" \
  --layout reference-walk \
  --reference-walk GRCh38 \
  --highlight-walk SAMPLE#1 \
  --highlight-walk SAMPLE#2 \
  --fade-opacity 0.15 \
  --title "${GENE}"
```

Manuscript schematic workflow after editing a layout:

```bash
node --check scripts/pangene_gfa_to_svg.js

scripts/pangene_gfa_to_svg.js \
  --input INPUT.gfa \
  --output manuscript/figs/vector/FIGURE.events.white.svg \
  --layout LAYOUT_NAME \
  --background white \
  --title "FIGURE TITLE"

xmllint --noout manuscript/figs/vector/FIGURE.events.white.svg

rsvg-convert -f pdf \
  manuscript/figs/vector/FIGURE.events.white.svg \
  -o manuscript/figs/vector/FIGURE.events.white.pdf

pdfinfo manuscript/figs/vector/FIGURE.events.white.pdf
```

For event schematics, check that `pdfinfo` reports a page width of `504 pt`, which is 7 inches.

### Prompt Templates

Use these prompts when asking an assistant to update the existing figures or create new ones.

Template: adjust an existing event figure

```text
Please update the existing Pangene event schematic in scripts/pangene_gfa_to_svg.js.

Target layout: [event-tracks or pms2-events]
Input GFA: [path/to/input.gfa]
Outputs to regenerate:
- manuscript/figs/vector/[NAME].events.svg
- manuscript/figs/vector/[NAME].events.white.svg
- matching PDFs

Requested visual changes:
- [e.g. move label X above track Y]
- [e.g. collapse region A-B]
- [e.g. change annotation text to "..."]

Keep the manuscript sizing constraints:
- SVG physical width no more than 7 inches
- visible text at least 10 pt
- no individual gene labels inside arrow blocks
- regenerate PDFs and validate SVG XML
```

Template: make a new reference-walk figure from a Pangene GFA

```text
Please generate a reference-walk SVG/PDF figure from this Pangene GFA:

Input GFA: [path/to/gene.gfa]
Gene/region name: [NAME]
Reference walk: [GRCh38, CHM13#0, or other]
Highlighted walks: [SAMPLE#1], [SAMPLE#2]
Main manuscript point: [one or two sentences]

Start with --layout reference-walk unless the manuscript point requires a custom schematic.
If a custom schematic is needed, explain why and then implement the smallest gene-specific layout.
Regenerate SVG and PDF outputs under manuscript/figs/vector/ and validate them.
```

Template: make a new custom event schematic

```text
Please create a new custom event schematic by extending scripts/pangene_gfa_to_svg.js.

Input GFA: [path/to/gene.gfa]
Output prefix: manuscript/figs/vector/[NAME]
Tracks/walks to show: [GRCh38, HG002#1, HG002#2, etc.]
Reference order: [walk name]
Events to emphasize:
- [event 1, including genes and size if known]
- [event 2]
- [benchmark/alignment interpretation to annotate]
Regions that can be collapsed:
- [left flank]
- [right flank]

Use GFA gene colors where possible. Verify expected walks and required genes are present.
Keep the figure manuscript-ready: width <= 7 inches, text >= 10 pt, no individual gene labels inside arrows, editable SVG shapes, and regenerated PDFs.
Update docs/pangene-vector-export.md with the exact command and any caveats.
```

Template: review or polish before manuscript use

```text
Please review the Pangene figure script and generated outputs for manuscript use.

Check:
- script syntax
- SVG XML validity
- PDF generation
- PDF physical width
- no unintended per-gene labels in event schematics
- text is readable and not overlapping
- documentation has reproducible commands

Summarize any remaining caveats and provide a short daily-note summary.
```

Convert to PDF:

```bash
inkscape manuscript/figs/vector/SULT1A1.pangene.svg \
  --export-filename=manuscript/figs/vector/SULT1A1.pangene.pdf
```

If Inkscape is not installed but `rsvg-convert` is available:

```bash
rsvg-convert -f pdf \
  manuscript/figs/vector/SULT1A1.pangene.svg \
  -o manuscript/figs/vector/SULT1A1.pangene.pdf
```

To apply this to another gene, change the `gene`, output basename, and optional `title`:

```bash
GENE=CYP2D6

curl -L \
  "https://pangene.bioinweb.org/view?graph=human472-1.1a2&gene=${GENE}&step=3&ori=" \
  -o "work/pangene/${GENE}.pangene.html"

scripts/extract_pangene_gfa.py \
  --input "work/pangene/${GENE}.pangene.html" \
  --output "work/pangene/${GENE}.gfa"

scripts/pangene_gfa_to_svg.js \
  --input "work/pangene/${GENE}.gfa" \
  --output "manuscript/figs/vector/${GENE}.pangene.svg" \
  --title "${GENE}"
```

For gene lists or names with special URL characters, URL-encode the `gene=` value before fetching the Pangene page.

Troubleshooting:

- If an SVG editor shows labels but not gene blocks, regenerate the SVG with the current script. Older script output used CSS `hsl(...)` colors that some editors did not render reliably.
- If the blocks are visible but too thin for a publication figure, use `--block-height 14` or another modest value and regenerate the PDF from that SVG.

Pros:

- Most reproducible for the paper.
- Does not require maintaining a forked web server.
- Produces real SVG with editable text and shapes.

Cons:

- Requires two small helper scripts.
- If upstream layout code changes, the script needs to be refreshed.

## Approach 3: Extract The Subgraph And Render With Existing Vector Tools

This avoids modifying Pangene/GFA viewer code. Extract the same local GFA subgraph, then use a vector-capable graph renderer. This is useful if exact Pangene styling is less important than a clean vector diagram.

Get the local GFA exactly as in Approach 2:

```bash
curl -L \
  "https://pangene.bioinweb.org/view?graph=human472-1.1a2&gene=SULT1A1&step=3&ori=" \
  -o work/pangene/SULT1A1.pangene.html

scripts/extract_pangene_gfa.py \
  --input work/pangene/SULT1A1.pangene.html \
  --output work/pangene/SULT1A1.gfa
```

Option A: open the GFA in BandageNG and export SVG/PDF from the GUI:

```bash
BandageNG work/pangene/SULT1A1.gfa
```

Option B: convert GFA links to DOT and render with Graphviz:

```bash
python3 scripts/gfa_to_dot.py \
  work/pangene/SULT1A1.gfa \
  > work/pangene/SULT1A1.dot

dot -Tsvg work/pangene/SULT1A1.dot \
  > manuscript/figs/vector/SULT1A1.graphviz.svg

dot -Tpdf work/pangene/SULT1A1.dot \
  > manuscript/figs/vector/SULT1A1.graphviz.pdf
```

Pros:

- Minimal Pangene-specific code.
- Quickly gives editable vector output.
- Good fallback for exploratory figures.

Cons:

- Layout and colors will differ from the Pangene webapp.
- Some biological visual conventions from Pangene, such as arrow-block styling and loss-of-function hollow blocks, need custom styling.

## Bottom Line

Yes, modifying the code to produce SVG/PDF is feasible. The graph is rasterized only at the final browser drawing step. The upstream data and layout are already vector-friendly.

For this manuscript, Approach 2 is implemented with:

- `scripts/extract_pangene_gfa.py`
- `scripts/pangene_gfa_to_svg.js`
- `manuscript/figs/vector/SULT1A1.pangene.svg`
- `manuscript/figs/vector/SULT1A1.pangene.white.svg`
- `manuscript/figs/vector/SULT1A1.pangene.pdf`
- `manuscript/figs/vector/SULT1A1.hg002-linear.svg`
- `manuscript/figs/vector/SULT1A1.hg002-linear.white.svg`
- `manuscript/figs/vector/SULT1A1.hg002-linear.pdf`
- `manuscript/figs/vector/SULT1A1.events.svg`
- `manuscript/figs/vector/SULT1A1.events.pdf`
- `manuscript/figs/vector/SULT1A1.events.white.svg`
- `manuscript/figs/vector/SULT1A1.events.white.pdf`
- `manuscript/figs/vector/PMS2.events.svg`
- `manuscript/figs/vector/PMS2.events.pdf`
- `manuscript/figs/vector/PMS2.events.white.svg`
- `manuscript/figs/vector/PMS2.events.white.pdf`

If the result is useful beyond this paper, promote the same SVG drawing functions into `gfatools/js/gfa-plot.js` as Approach 1.
