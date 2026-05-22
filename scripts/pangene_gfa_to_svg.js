#!/usr/bin/env node
"use strict";

/*
 * Render a Pangene/gfa-server local GFA subgraph to SVG.
 *
 * The parser and layout logic are adapted from gfatools js/gfa.js and
 * js/gfa-plot.js. This script replaces the final canvas drawing calls with
 * SVG elements so the result is editable in vector graphics tools.
 */

const fs = require("fs");
const path = require("path");

function usage() {
  console.error(`Usage:
  node scripts/pangene_gfa_to_svg.js \\
    --input work/pangene/SULT1A1.gfa \\
    --output manuscript/figs/vector/SULT1A1.pangene.svg

Options:
  --input PATH       Input local GFA subgraph
  --output PATH      Output SVG
  --layout STR       Layout: pangene, reference-walk, event-tracks, or pms2-events [pangene]
  --reference-walk STR
                    Reference walk for reference-walk/event layouts [GRCh38]
  --highlight-walk STR
                    Walk to highlight; repeatable. Defaults to HG002#1/HG002#2
                    for reference-walk and event-tracks. In pms2-events, an
                    inverted walk is selected automatically unless provided.
  --fade-opacity NUM Opacity for genes absent from highlighted walks [0.2]
  --label STR        Node label: name or length [name]
  --font-size INT    Label font size in px [9]
  --block-height INT Full gene block height in px [10]
  --background STR   Optional SVG background color, e.g. white or none [none]
  --title STR        Optional SVG title

Notes:
  pangene and reference-walk are general-purpose layouts.
  event-tracks and pms2-events are manuscript schematics capped at 7 inches
  wide, with visible text at 10 pt or larger and unlabeled gene blocks.
`);
}

function parseArgs(argv) {
  const args = {
    layout: "pangene",
    referenceWalk: null,
    highlightWalks: [],
    fadeOpacity: 0.2,
    label: "name",
    fontSize: 9,
    blockHeight: 10,
    background: "none",
    title: null,
  };

  for (let i = 2; i < argv.length; ++i) {
    const a = argv[i];
    const next = () => {
      if (i + 1 >= argv.length) {
        throw new Error(`Missing value for ${a}`);
      }
      return argv[++i];
    };

    if (a === "--input") args.input = next();
    else if (a === "--output") args.output = next();
    else if (a === "--layout") args.layout = next();
    else if (a === "--reference-walk") args.referenceWalk = next();
    else if (a === "--highlight-walk") args.highlightWalks.push(next());
    else if (a === "--fade-opacity") args.fadeOpacity = Number.parseFloat(next());
    else if (a === "--label") args.label = next();
    else if (a === "--font-size") args.fontSize = Number.parseInt(next(), 10);
    else if (a === "--block-height") args.blockHeight = Number.parseInt(next(), 10);
    else if (a === "--background") args.background = next();
    else if (a === "--title") args.title = next();
    else if (a === "-h" || a === "--help") {
      usage();
      process.exit(0);
    } else {
      throw new Error(`Unknown option: ${a}`);
    }
  }

  if (!args.input || !args.output) {
    usage();
    process.exit(2);
  }
  if (!["pangene", "reference-walk", "event-tracks", "pms2-events"].includes(args.layout)) {
    throw new Error("--layout must be 'pangene', 'reference-walk', 'event-tracks', or 'pms2-events'");
  }
  if (args.layout === "reference-walk" || args.layout === "event-tracks" || args.layout === "pms2-events") {
    if (args.referenceWalk == null) args.referenceWalk = "GRCh38";
  }
  if (args.layout === "reference-walk" || args.layout === "event-tracks") {
    if (args.highlightWalks.length === 0) args.highlightWalks = ["HG002#1", "HG002#2"];
  }
  if (!Number.isFinite(args.fadeOpacity) || args.fadeOpacity < 0 || args.fadeOpacity > 1) {
    throw new Error("--fade-opacity must be a number between 0 and 1");
  }
  if (!["name", "length"].includes(args.label)) {
    throw new Error("--label must be 'name' or 'length'");
  }
  if (!Number.isFinite(args.fontSize) || args.fontSize <= 0) {
    throw new Error("--font-size must be a positive integer");
  }
  if (!Number.isFinite(args.blockHeight) || args.blockHeight <= 0) {
    throw new Error("--block-height must be a positive integer");
  }
  return args;
}

function plotConf(args) {
  return {
    label: args.label,
    font_size: args.fontSize,
    min_len: 10,
    scale: 10,
    h_arrow: args.blockHeight / 2,
    xskip: 15,
    yskip: 30,
  };
}

function segAdd(g, name) {
  if (Object.prototype.hasOwnProperty.call(g.segname, name)) {
    return g.segname[name];
  }
  const sid = g.seg.length;
  g.segname[name] = sid;
  g.seg.push({ name, len: -1, sname: null, soff: -1, rank: -1, color: null });
  return sid;
}

function indexGraph(g) {
  const nVtx = g.seg.length * 2;
  for (let v = 0; v < nVtx; ++v) g.idx[v] = { o: 0, n: 0 };

  g.arc.sort((a, b) => a.v - b.v);
  let st = 0;
  for (let i = 1; i <= g.arc.length; ++i) {
    if (i === g.arc.length || g.arc[i].v !== g.arc[st].v) {
      g.idx[g.arc[st].v] = { o: st, n: i - st };
      st = i;
    }
  }

  for (let v = 0; v < nVtx; ++v) {
    const ov = g.idx[v].o;
    const nv = g.idx[v].n;
    let i0 = -1;
    let n0 = 0;
    for (let i = 0; i < nv; ++i) {
      if (g.arc[ov + i].rank === 0) {
        ++n0;
        i0 = i;
      }
    }
    if (n0 > 1) g.err |= 2;
    if (i0 > 0) {
      const tmp = g.arc[ov];
      g.arc[ov] = g.arc[ov + i0];
      g.arc[ov + i0] = tmp;
    }
  }
}

function hslToHex(h, s, l) {
  h = ((h % 360) + 360) % 360;
  s /= 100;
  l /= 100;

  const c = (1 - Math.abs(2 * l - 1)) * s;
  const hp = h / 60;
  const x = c * (1 - Math.abs((hp % 2) - 1));
  let r = 0;
  let g = 0;
  let b = 0;

  if (hp >= 0 && hp < 1) [r, g, b] = [c, x, 0];
  else if (hp < 2) [r, g, b] = [x, c, 0];
  else if (hp < 3) [r, g, b] = [0, c, x];
  else if (hp < 4) [r, g, b] = [0, x, c];
  else if (hp < 5) [r, g, b] = [x, 0, c];
  else [r, g, b] = [c, 0, x];

  const m = l - c / 2;
  const toHex = (v) => {
    const n = Math.round((v + m) * 255);
    return n.toString(16).padStart(2, "0");
  };
  return `#${toHex(r)}${toHex(g)}${toHex(b)}`;
}

function assignColor(g) {
  for (let i = 0; i < g.seg.length; ++i) {
    g.seg[i].color = hslToHex(i * 137.508, 65, 50);
  }
}

function parseGfa(str) {
  const g = { seg: [], arc: [], segname: {}, idx: [], walk: [], err: 0 };
  const lines = str.split(/\r?\n/);
  const reCigar = /(\d+)([MIDSN])/g;
  const reWalk = /([><])([^\s><]+)/g;

  for (const line of lines) {
    if (line.length < 5) continue;
    const t = line.split("\t");
    let m;

    if (t[0] === "S") {
      const sid = segAdd(g, t[1]);
      const s = g.seg[sid];
      if (t[2] !== "*") s.len = t[2].length;
      for (let j = 3; j < t.length; ++j) {
        m = /^(LN:i|SN:Z|SO:i|SR:i):(\S+)/.exec(t[j]);
        if (m == null) continue;
        if (m[1] === "LN:i") s.len = Number.parseInt(m[2], 10);
        else if (m[1] === "SN:Z") s.sname = m[2];
        else if (m[1] === "SO:i") s.soff = Number.parseInt(m[2], 10);
        else if (m[1] === "SR:i") s.rank = Number.parseInt(m[2], 10);
      }
    } else if (t[0] === "L") {
      if (t.length < 5 || !["+", "-"].includes(t[2]) || !["+", "-"].includes(t[4])) {
        continue;
      }
      const sid1 = segAdd(g, t[1]);
      const sid2 = segAdd(g, t[3]);
      const v = (sid1 << 1) | (t[2] === "+" ? 0 : 1);
      const w = (sid2 << 1) | (t[4] === "+" ? 0 : 1);
      let ov = 0;
      let ow = 0;
      let rank = -1;

      for (let j = 6; j < t.length; ++j) {
        m = /^(SR:i):(\S+)/.exec(t[j]);
        if (m != null) rank = Number.parseInt(m[2], 10);
      }
      if (t.length >= 6) {
        reCigar.lastIndex = 0;
        while ((m = reCigar.exec(t[5])) != null) {
          const n = Number.parseInt(m[1], 10);
          if (m[2] === "M" || m[2] === "D" || m[2] === "N") ov += n;
          if (m[2] === "M" || m[2] === "I" || m[2] === "S") ow += n;
        }
      }
      g.arc.push({ v, w, ov, ow, rank, ori: true });
      g.arc.push({ v: w ^ 1, w: v ^ 1, ov: ow, ow: ov, rank, ori: false });
    } else if (t[0] === "W") {
      if (t.length < 7) continue;
      const walk = {
        asm: `${t[1]}#${t[2]}`,
        sample: t[1],
        hap: Number.parseInt(t[2], 10),
        sname: t[3],
        lineName: t[3],
        st: -1,
        en: -1,
        v: [],
        lof: [],
      };
      if (t[4] !== "*") walk.st = Number.parseInt(t[4], 10);
      if (t[5] !== "*") walk.en = Number.parseInt(t[5], 10);
      reWalk.lastIndex = 0;
      while ((m = reWalk.exec(t[6])) != null) {
        if (Object.prototype.hasOwnProperty.call(g.segname, m[2])) {
          const sid = g.segname[m[2]];
          walk.v.push((sid << 1) | (m[1] === ">" ? 0 : 1));
        }
      }
      g.walk.push(walk);
    }
  }

  for (let i = 0; i < g.seg.length; ++i) {
    if (g.seg[i].len < 0) g.err |= 1;
  }
  if (g.seg.length === 0) {
    throw new Error("input GFA has no segment records");
  }
  indexGraph(g);
  assignColor(g);
  return g;
}

function scc1Aux(g) {
  const nVtx = g.seg.length * 2;
  const aux = { a: [], index: 0 };
  for (let i = 0; i < nVtx; ++i) {
    aux.a.push({ index: -1, start: -1, low: 0, i: -1, stack: false });
  }
  return aux;
}

function scc1(g, aux, v0) {
  const sub = { v: [], a: [] };
  const ds = [[v0, 0]];
  const ts = [];

  while (ds.length > 0) {
    let [v, i] = ds.pop();
    if (i === 0) {
      aux.a[v].low = aux.a[v].index = aux.index++;
      aux.a[v].stack = true;
      ts.push(v);
    }
    const nv = g.idx[v].n;
    if (i === nv) {
      if (aux.a[v].low === aux.a[v].index) {
        while (ts.length > 0) {
          const w = ts.pop();
          sub.v.push({ v: w, off: 0, n: 0 });
          aux.a[w].stack = false;
          if (w === v) break;
        }
      }
      if (ds.length > 0) {
        const w = v;
        v = ds[ds.length - 1][0];
        aux.a[v].low = Math.min(aux.a[v].low, aux.a[w].low);
      }
    } else {
      const w = g.arc[g.idx[v].o + i].w;
      ds.push([v, i + 1]);
      if (aux.a[w].index === -1 && aux.a[w ^ 1].stack === false) {
        ds.push([w, 0]);
      } else if (aux.a[w].stack) {
        aux.a[v].low = Math.min(aux.a[v].low, aux.a[w].index);
      }
    }
  }

  sub.v.reverse();
  for (let k = 0; k < sub.v.length; ++k) {
    aux.a[sub.v[k].v].start = v0;
    aux.a[sub.v[k].v].i = k;
  }
  for (let k = 0; k < sub.v.length; ++k) {
    const o0 = sub.a.length;
    const v = sub.v[k].v;
    const nv = g.idx[v].n;
    const ov = g.idx[v].o;
    for (let i = 0; i < nv; ++i) {
      const a = g.arc[ov + i];
      if (aux.a[a.w].start === v0) {
        sub.a.push({ i: aux.a[a.w].i, arc_off: ov + i, rank: a.rank });
      }
    }
    sub.v[k].off = o0;
    sub.v[k].n = sub.a.length - o0;
    if (sub.v[k].n > 1) {
      const sorted = sub.a.slice(o0).sort((x, y) => x.i - y.i);
      for (let i = 0; i < sorted.length; ++i) sub.a[o0 + i] = sorted[i];
    }
  }
  return sub;
}

function rankToColor(g) {
  const colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#a65628", "#f781bf"];
  const ranks = [];
  for (const s of g.seg) if (s.rank >= 0) ranks.push(s.rank);
  for (const a of g.arc) if (a.rank >= 0) ranks.push(a.rank);
  if (ranks.length === 0) return {};

  ranks.sort((a, b) => a - b);
  const uniq = [ranks[0]];
  for (let i = 1; i < ranks.length; ++i) {
    if (ranks[i] !== ranks[i - 1]) uniq.push(ranks[i]);
  }
  const r2c = {};
  for (let i = 0; i < uniq.length && i < colors.length; ++i) {
    r2c[uniq[i]] = colors[i];
  }
  return r2c;
}

function findV0(g) {
  const nVtx = g.seg.length * 2;
  let v0Ref = -1;
  let v0Src = -1;

  for (let v = 0; v < nVtx; ++v) {
    const s = g.seg[v >> 1];
    if (s.rank === 0) {
      if (s.snid < 0 || s.soff < 0) continue;
      if (v0Ref < 0 || s.soff < g.seg[v0Ref >> 1].soff) v0Ref = v;
    }
    if (g.idx[v ^ 1].n === 0 && v0Src < 0) v0Src = v;
  }
  return v0Ref >= 0 ? v0Ref : v0Src >= 0 ? v0Src : 0;
}

function calLength(len, minLen, scale) {
  return Math.floor(minLen + (Math.log(len + 1) / Math.log(10)) * scale + 0.499);
}

function calPos(conf, g) {
  const v0 = findV0(g);
  const aux = scc1Aux(g);
  const sub = scc1(g, aux, v0);
  const pred = [];
  for (let i = 0; i < sub.v.length; ++i) pred[i] = [];
  for (let i = 0; i < sub.v.length; ++i) {
    for (let j = 0; j < sub.v[i].n; ++j) {
      pred[sub.a[sub.v[i].off + j].i].push(i);
    }
  }

  const levelMax = [];
  const pos = [];
  for (let i = 0; i < sub.v.length; ++i) {
    pos[i] = { level: -1, start: -1, len: 0 };
    pos[i].len = calLength(g.seg[sub.v[i].v >> 1].len, conf.min_len, conf.scale);
  }

  for (let i = 0; i < sub.v.length; ++i) {
    if (pred[i].length === 0) {
      pos[i].start = 0;
      pos[i].level = levelMax.length;
    } else {
      let maxEnd = -1;
      const pl = [];
      for (let j = 0; j < levelMax.length; ++j) pl[j] = { cnt: 0, i: -1, end: -1 };

      for (let j = 0; j < pred[i].length; ++j) {
        const p = pos[pred[i][j]];
        if (p.level >= 0) {
          const end = p.start + p.len;
          if (end > maxEnd) maxEnd = end;
          pl[p.level].i = pred[i][j];
          pl[p.level].end = end;
          ++pl[p.level].cnt;
        }
      }

      pos[i].start = maxEnd + conf.xskip;
      let l;
      for (l = 0; l < levelMax.length; ++l) {
        if (pl[l].cnt > 1 || levelMax[l] + conf.xskip > pos[i].start) continue;
        if (pl[l].cnt === 0) break;
        if (pl[l].end === levelMax[l]) break;
      }
      pos[i].level = l;
    }
    levelMax[pos[i].level] = pos[i].start + pos[i].len;
  }

  return [sub, pos];
}

function xmlEscape(s) {
  return String(s)
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;");
}

function physicalSizeAttrs(width, height, maxWidthIn = 7) {
  const heightIn = (maxWidthIn * height) / width;
  const h = heightIn.toFixed(3).replace(/0+$/, "").replace(/\.$/, "");
  return `width="${maxWidthIn}in" height="${h}in"`;
}

function arrowPoints(x, y, len, w, rev) {
  if (rev) {
    return [
      [x, y],
      [x + w, y - w],
      [x + w + len, y - w],
      [x + len, y],
      [x + w + len, y + w],
      [x + w, y + w],
    ];
  }
  return [
    [x, y],
    [x - w, y - w],
    [x - w + len, y - w],
    [x + len, y],
    [x - w + len, y + w],
    [x - w, y + w],
  ];
}

function svgArrow(out, x, y, len, w, rev, text, fs, stroke, fill) {
  const pts = arrowPoints(x, y, len, w, rev).map((p) => p.join(",")).join(" ");
  const shapeStroke = stroke || fill || "#000000";
  const attrs = [
    `points="${pts}"`,
    `fill="${xmlEscape(fill || "none")}"`,
    `stroke="${xmlEscape(shapeStroke)}"`,
    `stroke-width="1"`,
    `vector-effect="non-scaling-stroke"`,
  ];
  out.push(`  <polygon ${attrs.join(" ")} />`);

  if (text != null) {
    const textColor = stroke || fill || "#000000";
    out.push(
      `  <text x="${x + len / 2}" y="${y - w - 2}" text-anchor="middle" ` +
        `font-family="monospace" font-size="${fs}" fill="${xmlEscape(textColor)}">` +
        `${xmlEscape(text)}</text>`
    );
  }
}

function walkDisplayName(w) {
  return `${w.sample}#${w.hap}`;
}

function findWalk(g, label) {
  const exact = g.walk.filter(
    (w) => walkDisplayName(w) === label || w.asm === label || w.lineName === label
  );
  if (exact.length === 1) return exact[0];

  if (exact.length > 1) {
    throw new Error(
      `walk label '${label}' matched multiple walks: ${exact.map(walkDisplayName).join(", ")}`
    );
  }

  const bySample = g.walk.filter((w) => w.sample === label);
  if (bySample.length === 1) return bySample[0];
  if (bySample.length > 1) {
    throw new Error(
      `walk label '${label}' matched multiple haplotypes; use one of: ` +
        bySample.map(walkDisplayName).join(", ")
    );
  }

  throw new Error(`failed to find walk '${label}' in GFA`);
}

function uniqueWalkGenes(walks) {
  const genes = new Set();
  for (const walk of walks) {
    for (const v of walk.v) genes.add(v >> 1);
  }
  return genes;
}

function itemCenter(item, yOffset = 0) {
  return {
    x: item.cx_st + item.len / 2,
    y: item.cy + yOffset,
  };
}

function slugId(s) {
  return String(s).replace(/[^A-Za-z0-9_.-]+/g, "_");
}

function makeLayoutItem({ key, sid, v, x, y, len, lane, occurrence }) {
  return {
    key,
    sid,
    v,
    len,
    lane,
    occurrence,
    cx_st: x,
    cx_en: x + len,
    cy: y,
  };
}

function placeRemainingSegments(g, conf, placedItems, y, maxWidth) {
  const placedSids = new Set(placedItems.map((item) => item.sid));
  const items = [];
  let x = conf.xskip;

  for (let sid = 0; sid < g.seg.length; ++sid) {
    if (placedSids.has(sid)) continue;
    const len = calLength(g.seg[sid].len, conf.min_len, conf.scale);
    if (x > conf.xskip && x + len + conf.xskip > maxWidth) {
      x = conf.xskip;
      y += conf.yskip;
    }
    items.push(
      makeLayoutItem({
        key: `unhighlighted:${sid}`,
        sid,
        v: sid << 1,
        x,
        y,
        len,
        lane: "unhighlighted",
        occurrence: 0,
      })
    );
    x += len + conf.xskip;
  }

  return items;
}

function layoutReferenceGenes(g, referenceWalk, conf) {
  const items = [];
  const refByIndex = [];
  const refIndicesBySid = new Map();
  let x = conf.xskip;
  const y = 90;

  for (let i = 0; i < referenceWalk.v.length; ++i) {
    const v = referenceWalk.v[i];
    const sid = v >> 1;
    const len = calLength(g.seg[sid].len, conf.min_len, conf.scale);
    const item = makeLayoutItem({
      key: `ref:${i}`,
      sid,
      v,
      x,
      y,
      len,
      lane: "reference",
      occurrence: i,
    });
    items.push(item);
    refByIndex[i] = item;
    if (!refIndicesBySid.has(sid)) refIndicesBySid.set(sid, []);
    refIndicesBySid.get(sid).push(i);
    x += len + conf.xskip;
  }

  return { items, refByIndex, refIndicesBySid };
}

function mapWalkToReference(walk, refIndicesBySid) {
  const mapped = [];
  let cursor = 0;

  for (let i = 0; i < walk.v.length; ++i) {
    const sid = walk.v[i] >> 1;
    const candidates = refIndicesBySid.get(sid) || [];
    let refIndex = null;
    for (const idx of candidates) {
      if (idx >= cursor) {
        refIndex = idx;
        break;
      }
    }
    if (refIndex != null) cursor = refIndex + 1;
    mapped.push({ walkIndex: i, v: walk.v[i], sid, refIndex, itemKey: null });
  }

  return mapped;
}

function placeUnmatchedRun(g, conf, run, prevRef, nextRef, laneY, laneName, keyPrefix) {
  const totalLen =
    run.reduce((sum, step) => sum + calLength(g.seg[step.sid].len, conf.min_len, conf.scale), 0) +
    Math.max(0, run.length - 1) * conf.xskip;
  const left = prevRef ? prevRef.cx_en + conf.xskip : conf.xskip;
  const right = nextRef ? nextRef.cx_st - conf.xskip : left + totalLen;
  let x = left;

  if (right > left && right - left > totalLen) {
    x = left + (right - left - totalLen) / 2;
  }

  return run.map((step, i) => {
    const len = calLength(g.seg[step.sid].len, conf.min_len, conf.scale);
    const item = makeLayoutItem({
      key: `${keyPrefix}:alt:${step.walkIndex}`,
      sid: step.sid,
      v: step.v,
      x,
      y: laneY,
      len,
      lane: laneName,
      occurrence: step.walkIndex,
    });
    x += len + conf.xskip;
    return item;
  });
}

function placeHighlightedWalk(g, conf, walk, refByIndex, refIndicesBySid, laneY, laneName, keyPrefix) {
  const mapped = mapWalkToReference(walk, refIndicesBySid);
  const items = [];
  const stepItems = [];

  for (let i = 0; i < mapped.length; ) {
    if (mapped[i].refIndex != null) {
      const item = refByIndex[mapped[i].refIndex];
      mapped[i].itemKey = item.key;
      stepItems[i] = item;
      ++i;
      continue;
    }

    const start = i;
    while (i < mapped.length && mapped[i].refIndex == null) ++i;
    const run = mapped.slice(start, i);
    const prevMapped = mapped.slice(0, start).reverse().find((step) => step.refIndex != null);
    const nextMapped = mapped.slice(i).find((step) => step.refIndex != null);
    const prevRef = prevMapped ? refByIndex[prevMapped.refIndex] : null;
    const nextRef = nextMapped ? refByIndex[nextMapped.refIndex] : null;
    const runItems = placeUnmatchedRun(g, conf, run, prevRef, nextRef, laneY, laneName, keyPrefix);
    for (let j = 0; j < runItems.length; ++j) {
      items.push(runItems[j]);
      mapped[start + j].itemKey = runItems[j].key;
      stepItems[start + j] = runItems[j];
    }
  }

  return { mapped, items, stepItems };
}

function renderPolyline(out, id, points, color, width, opacity = 1) {
  if (points.length < 2) return;
  out.push(
    `  <polyline id="${xmlEscape(slugId(id))}" points="${points.map((p) => `${p.x},${p.y}`).join(" ")}" ` +
      `fill="none" stroke="${xmlEscape(color)}" stroke-width="${width}" ` +
      `stroke-linecap="round" stroke-linejoin="round" opacity="${opacity}" ` +
      `vector-effect="non-scaling-stroke" />`
  );
}

function renderGeneItem(out, g, conf, item, highlightedGenes, fadeOpacity) {
  const s = g.seg[item.sid];
  const label = conf.label === "name" ? s.name : s.len;
  const opacity = highlightedGenes.has(item.sid) ? 1 : fadeOpacity;
  out.push(
    `  <g class="gene" data-gene="${xmlEscape(s.name)}" data-lane="${xmlEscape(item.lane)}" ` +
      `data-occurrence="${item.occurrence}" opacity="${opacity}">`
  );
  svgArrow(
    out,
    item.cx_st,
    item.cy,
    item.len,
    conf.h_arrow,
    Boolean(item.v & 1),
    label,
    conf.font_size,
    null,
    s.color
  );
  out.push(`  </g>`);
}

function renderReferenceWalkSvg(conf, g, options) {
  const referenceWalk = findWalk(g, options.referenceWalk);
  const highlightWalks = options.highlightWalks.map((label) => findWalk(g, label));
  const highlightedGenes = uniqueWalkGenes(highlightWalks);
  const refLayout = layoutReferenceGenes(g, referenceWalk, conf);
  const allItems = [...refLayout.items];
  const pathColors = ["#0072B2", "#D55E00", "#009E73", "#CC79A7"];
  const pathOffsets = [-6, 6, -12, 12];
  const highlightedPaths = [];

  for (let i = 0; i < highlightWalks.length; ++i) {
    const walk = highlightWalks[i];
    const laneY = i === 0 ? 45 : i === 1 ? 135 : 45 + i * 30;
    const laneName = `highlight_${i + 1}`;
    const placed = placeHighlightedWalk(
      g,
      conf,
      walk,
      refLayout.refByIndex,
      refLayout.refIndicesBySid,
      laneY,
      laneName,
      `highlight:${i}`
    );
    allItems.push(...placed.items);
    highlightedPaths.push({
      label: walkDisplayName(walk),
      color: pathColors[i % pathColors.length],
      points: placed.stepItems.filter(Boolean).map((item) => itemCenter(item, pathOffsets[i] || 0)),
    });
  }

  const lastItem = allItems.length > 0 ? allItems[allItems.length - 1] : null;
  allItems.push(
    ...placeRemainingSegments(g, conf, allItems, 180, Math.max(600, lastItem ? lastItem.cx_en : 600))
  );

  let maxW = 0;
  let maxH = 0;
  for (const item of allItems) {
    maxW = Math.max(maxW, item.cx_en + conf.xskip);
    maxH = Math.max(maxH, item.cy + conf.yskip);
  }

  const out = [];
  out.push(`<?xml version="1.0" encoding="UTF-8"?>`);
  out.push(
    `<svg xmlns="http://www.w3.org/2000/svg" width="${maxW}" height="${maxH}" ` +
      `viewBox="0 0 ${maxW} ${maxH}">`
  );
  if (options.title) out.push(`  <title>${xmlEscape(options.title)}</title>`);
  if (options.background !== "none") {
    out.push(`  <rect width="100%" height="100%" fill="${xmlEscape(options.background)}" />`);
  }

  const refPoints = refLayout.items.map((item) => itemCenter(item, 0));
  out.push(`  <g id="reference-paths">`);
  renderPolyline(out, `walk-${walkDisplayName(referenceWalk)}`, refPoints, "#606060", 1.5, 0.45);
  out.push(`  </g>`);

  out.push(`  <g id="highlight-paths">`);
  for (const pathInfo of highlightedPaths) {
    renderPolyline(out, `walk-${pathInfo.label}`, pathInfo.points, pathInfo.color, 3, 0.95);
  }
  out.push(`  </g>`);

  out.push(`  <g id="genes">`);
  const sortedItems = allItems.slice().sort((a, b) => a.cy - b.cy || a.cx_st - b.cx_st);
  for (const item of sortedItems) {
    renderGeneItem(out, g, conf, item, highlightedGenes, options.fadeOpacity);
  }
  out.push(`  </g>`);
  out.push(`</svg>`);
  return out.join("\n") + "\n";
}

function geneSid(g, name) {
  if (!Object.prototype.hasOwnProperty.call(g.segname, name)) {
    throw new Error(`event-tracks layout requires gene '${name}' in the input GFA`);
  }
  return g.segname[name];
}

function geneColor(g, name) {
  return g.seg[geneSid(g, name)].color || "#808080";
}

function addText(out, x, y, text, attrs = "") {
  out.push(`  <text x="${x}" y="${y}" ${attrs}>${xmlEscape(text)}</text>`);
}

function addLine(out, x1, y1, x2, y2, attrs = "") {
  out.push(`  <line x1="${x1}" y1="${y1}" x2="${x2}" y2="${y2}" ${attrs} />`);
}

function addCollapsedBlock(out, x, y, width, label) {
  out.push(
    `  <g class="collapsed-block" data-label="${xmlEscape(label)}">` +
      `<rect x="${x}" y="${y - 10}" width="${width}" height="20" rx="3" ` +
      `fill="#f2f2f2" stroke="#9a9a9a" stroke-width="1" vector-effect="non-scaling-stroke" />` +
      `<text x="${x + width / 2}" y="${y - 16}" text-anchor="middle" font-family="monospace" ` +
      `font-size="10pt" fill="#606060">${xmlEscape(label)}</text>` +
      `</g>`
  );
}

function addGap(out, x, y, width, label) {
  out.push(
    `  <g class="event deletion" data-label="${xmlEscape(label)}">` +
      `<rect x="${x}" y="${y - 10}" width="${width}" height="20" fill="none" ` +
      `stroke="#4d4d4d" stroke-width="1.2" stroke-dasharray="4 3" ` +
      `vector-effect="non-scaling-stroke" />` +
      `<text x="${x + width / 2}" y="${y + 26}" text-anchor="middle" font-family="monospace" ` +
      `font-size="10pt" fill="#4d4d4d">${xmlEscape(label)}</text>` +
      `</g>`
  );
}

function addBracket(out, x1, x2, y, label, color = "#333333") {
  addLine(out, x1, y, x2, y, `stroke="${color}" stroke-width="1.2" vector-effect="non-scaling-stroke"`);
  addLine(out, x1, y, x1, y + 7, `stroke="${color}" stroke-width="1.2" vector-effect="non-scaling-stroke"`);
  addLine(out, x2, y, x2, y + 7, `stroke="${color}" stroke-width="1.2" vector-effect="non-scaling-stroke"`);
  addText(
    out,
    (x1 + x2) / 2,
    y - 5,
    label,
    `text-anchor="middle" font-family="monospace" font-size="10pt" fill="${color}"`
  );
}

function addEventGene(out, g, conf, name, x, y, width, rev = false) {
  const sid = geneSid(g, name);
  out.push(`  <g class="gene event-gene" data-gene="${xmlEscape(name)}">`);
  svgArrow(out, x, y, width, conf.h_arrow, rev, null, conf.font_size, null, g.seg[sid].color);
  out.push(`  </g>`);
}

function addCassette(out, g, conf, x, y, scale = 1) {
  const widths = {
    BOLA2B: 34 * scale,
    SLX1A: 30 * scale,
    SULT1A3: 42 * scale,
  };
  const gap = 4 * scale;
  addEventGene(out, g, conf, "BOLA2B", x, y, widths.BOLA2B, false);
  x += widths.BOLA2B + gap;
  addEventGene(out, g, conf, "SLX1A", x, y, widths.SLX1A, true);
  x += widths.SLX1A + gap;
  addEventGene(out, g, conf, "SULT1A3", x, y, widths.SULT1A3, true);
  return x + widths.SULT1A3;
}

function addCassetteSeries(out, g, conf, x, y, copies, extraStart = null) {
  const startX = x;
  const copyWidth = 34 + 4 + 30 + 4 + 42;
  for (let i = 0; i < copies; ++i) {
    addCassette(out, g, conf, x, y, 1);
    if (extraStart != null && i + 1 >= extraStart) {
      out.push(
        `  <rect x="${x - 3}" y="${y - 13}" width="${copyWidth + 6}" height="26" ` +
          `fill="none" stroke="#D55E00" stroke-width="1" stroke-dasharray="3 2" ` +
          `vector-effect="non-scaling-stroke" />`
      );
    }
    x += copyWidth + 10;
  }
  return { startX, endX: x - 10 };
}

function renderEventTracksSvg(conf, g, options) {
  // Resolve walks to fail early if the expected comparison tracks are absent.
  findWalk(g, options.referenceWalk);
  options.highlightWalks.forEach((label) => findWalk(g, label));

  const eventConf = { ...conf, h_arrow: Math.max(conf.h_arrow, 8) };
  const width = 1510;
  const height = 350;
  const out = [];
  const trackY = {
    GRCh38: 95,
    "HG002#1": 180,
    "HG002#2": 265,
  };
  const labelX = 20;
  const x = {
    left: 130,
    deletion: 270,
    middle: 440,
    sult: 610,
    right: 750,
    cassette: 900,
  };

  out.push(`<?xml version="1.0" encoding="UTF-8"?>`);
  out.push(`<svg xmlns="http://www.w3.org/2000/svg" ${physicalSizeAttrs(width, height)} viewBox="0 0 ${width} ${height}">`);
  if (options.title) out.push(`  <title>${xmlEscape(options.title)}</title>`);
  if (options.background !== "none") {
    out.push(`  <rect width="100%" height="100%" fill="${xmlEscape(options.background)}" />`);
  }

  addText(
    out,
    20,
    24,
    "SULT1A1 region: HG002 haplotype-specific events",
    `font-family="Helvetica,Arial,sans-serif" font-size="12pt" font-weight="700" fill="#222222"`
  );

  for (const [label, y] of Object.entries(trackY)) {
    addText(out, labelX, y + 4, label, `font-family="monospace" font-size="10pt" fill="#222222"`);
    addLine(out, 110, y, 1265, y, `stroke="#d0d0d0" stroke-width="1" vector-effect="non-scaling-stroke"`);
  }

  for (const y of Object.values(trackY)) {
    addCollapsedBlock(out, x.left, y, 100, "XPO6-SBK1");
    addCollapsedBlock(out, x.middle, y, 120, "NPIPB7-SULT1A2");
    addCollapsedBlock(out, x.right, y, 120, "NPIPB8-NPIPB11");
  }

  addEventGene(out, g, eventConf, "NPIPB6", x.deletion, trackY.GRCh38, 46, false);
  addEventGene(out, g, eventConf, "EIF3C", x.deletion + 58, trackY.GRCh38, 46, false);
  addEventGene(out, g, eventConf, "NPIPB6", x.deletion, trackY["HG002#1"], 46, false);
  addEventGene(out, g, eventConf, "EIF3C", x.deletion + 58, trackY["HG002#1"], 46, false);
  addGap(out, x.deletion, trackY["HG002#2"], 104, "114 kb deletion");
  addBracket(out, x.deletion, x.deletion + 104, 52, "maternal deletion: NPIPB6 + EIF3C", "#4d4d4d");

  addEventGene(out, g, eventConf, "SULT1A1", x.sult, trackY.GRCh38, 50, false);
  addEventGene(out, g, eventConf, "SULT1A1", x.sult - 28, trackY["HG002#1"], 50, true);
  addEventGene(out, g, eventConf, "SULT1A1", x.sult + 28, trackY["HG002#1"], 50, false);
  addEventGene(out, g, eventConf, "SULT1A1", x.sult, trackY["HG002#2"], 50, false);
  addBracket(out, x.sult - 28, x.sult + 78, 137, "paternal +SULT1A1, 13 kbp", geneColor(g, "SULT1A1"));

  addCassetteSeries(out, g, eventConf, x.cassette, trackY.GRCh38, 2);
  const h1Cass = addCassetteSeries(out, g, eventConf, x.cassette, trackY["HG002#1"], 4, 3);
  const h2Cass = addCassetteSeries(out, g, eventConf, x.cassette, trackY["HG002#2"], 4, 3);
  addBracket(out, h1Cass.startX + 2 * 124, h1Cass.endX, 137, "+2 tandem BOLA2B-SLX1A-SULT1A3 copies", "#D55E00");
  addBracket(out, h2Cass.startX + 2 * 124, h2Cass.endX, 222, "~102 kbp insertion in v5.0q SV benchmark", "#D55E00");

  out.push(
    `  <g class="event chm13-inversion" opacity="0.9">` +
      `<line x1="118" y1="325" x2="1425" y2="325" stroke="#6A51A3" stroke-width="2" ` +
      `stroke-dasharray="6 4" vector-effect="non-scaling-stroke" />` +
      `<text x="690" y="318" text-anchor="middle" font-family="monospace" font-size="10pt" ` +
      `fill="#6A51A3">common inversion relative to CHM13; excluded from CHM13-based benchmark</text>` +
      `</g>`
  );

  out.push(`</svg>`);
  return out.join("\n") + "\n";
}

function walkGeneEntries(g, walk) {
  return walk.v.map((v) => ({ name: g.seg[v >> 1].name, rev: (v & 1) !== 0 }));
}

function containsGeneOrder(entries, names) {
  for (let i = 0; i <= entries.length - names.length; ++i) {
    let ok = true;
    for (let j = 0; j < names.length; ++j) {
      if (entries[i + j].name !== names[j]) {
        ok = false;
        break;
      }
    }
    if (ok) return true;
  }
  return false;
}

function findPms2InvertedWalks(g, invertedBlock) {
  const out = [];
  for (const walk of g.walk) {
    if (walk.sample === "GRCh38" || walk.sample === "CHM13") continue;
    if (containsGeneOrder(walkGeneEntries(g, walk), invertedBlock)) out.push(walk);
  }
  return out;
}

function addPath(out, d, attrs = "") {
  out.push(`  <path d="${d}" ${attrs} />`);
}

function addPms2Gene(out, g, conf, entry, x, y, width, highlight = false) {
  const sid = geneSid(g, entry.name);
  const fill = g.seg[sid].color;
  const stroke = highlight ? "#111111" : null;
  out.push(`  <g class="gene pms2-event-gene" data-gene="${xmlEscape(entry.name)}">`);
  svgArrow(out, x, y, width, conf.h_arrow, entry.rev, null, conf.font_size, stroke, fill);
  out.push(`  </g>`);
}

function drawPms2GeneSeries(out, g, conf, entries, x, y, width, gap, highlightGenes = new Set()) {
  const positions = new Map();
  for (let i = 0; i < entries.length; ++i) {
    const itemX = x + i * (width + gap);
    addPms2Gene(out, g, conf, entries[i], itemX, y, width, highlightGenes.has(entries[i].name));
    positions.set(entries[i].name, { x: itemX, y, width });
  }
  return positions;
}

function findEntriesForGenes(g, walk, names) {
  const entries = walkGeneEntries(g, walk);
  const wanted = new Set(names);
  return entries.filter((entry) => wanted.has(entry.name));
}

function renderPms2EventsSvg(conf, g, options) {
  const eventConf = { ...conf, h_arrow: Math.max(conf.h_arrow, 8) };
  const referenceWalk = findWalk(g, options.referenceWalk || "GRCh38");
  const normalBlock = [
    "PMS2",
    "AIMP2",
    "ANKRD61",
    "EIF2AK1",
    "USP42",
    "CYTH3",
    "FAM220A",
    "SAGSIN1",
    "RAC1",
    "DAGLB",
    "KDELR2",
    "GRID2IP",
    "ZDHHC4",
    "INTS15",
    "ZNF853",
    "ZNF316",
    "ZNF12",
  ];
  const invertedBlock = normalBlock.slice().reverse();
  normalBlock.forEach((name) => geneSid(g, name));

  const invertedWalks = findPms2InvertedWalks(g, invertedBlock);
  if (invertedWalks.length === 0) {
    throw new Error("pms2-events layout did not find any walks with the expected PMS2 inversion");
  }

  let invertedWalk = null;
  if (options.highlightWalks.length > 0) {
    for (const label of options.highlightWalks) {
      const candidate = findWalk(g, label);
      if (containsGeneOrder(walkGeneEntries(g, candidate), invertedBlock)) {
        invertedWalk = candidate;
        break;
      }
    }
  }
  if (invertedWalk == null) invertedWalk = invertedWalks[0];

  const refEntries = findEntriesForGenes(g, referenceWalk, normalBlock);
  const invEntries = findEntriesForGenes(g, invertedWalk, normalBlock);
  if (refEntries.length !== normalBlock.length || invEntries.length !== normalBlock.length) {
    throw new Error("pms2-events layout could not extract the expected 17 PMS2-region genes");
  }

  const width = 1510;
  const height = 360;
  const out = [];
  const trackY = {
    ref: 105,
    inv: 225,
  };
  const labelX = 20;
  const leftX = 130;
  const blockX = 290;
  const geneWidth = 47;
  const gap = 4;
  const blockEnd = blockX + normalBlock.length * geneWidth + (normalBlock.length - 1) * gap;
  const rightFlankX = blockEnd + 95;
  const highlightGenes = new Set(["PMS2", "RSPH10B", "RSPH10B2", "CCZ1B"]);
  const invertedCount = invertedWalks.length;

  out.push(`<?xml version="1.0" encoding="UTF-8"?>`);
  out.push(`<svg xmlns="http://www.w3.org/2000/svg" ${physicalSizeAttrs(width, height)} viewBox="0 0 ${width} ${height}">`);
  if (options.title) out.push(`  <title>${xmlEscape(options.title)}</title>`);
  if (options.background !== "none") {
    out.push(`  <rect width="100%" height="100%" fill="${xmlEscape(options.background)}" />`);
  }

  addText(
    out,
    20,
    24,
    "PMS2 segmental-duplication-mediated inversion",
    `font-family="Helvetica,Arial,sans-serif" font-size="12pt" font-weight="700" fill="#222222"`
  );
  addText(
    out,
    20,
    43,
    `${invertedCount} of 472 assembled haplotypes carry the 17-gene inversion in human472-1.1a2`,
    `font-family="Helvetica,Arial,sans-serif" font-size="10pt" fill="#444444"`
  );

  addText(out, labelX, trackY.ref + 4, walkDisplayName(referenceWalk), `font-family="monospace" font-size="10pt" fill="#222222"`);
  addText(
    out,
    labelX,
    trackY.inv + 4,
    `${walkDisplayName(invertedWalk)} inversion`,
    `font-family="monospace" font-size="10pt" fill="#222222"`
  );
  addLine(out, 130, trackY.ref, 1430, trackY.ref, `stroke="#d0d0d0" stroke-width="1" vector-effect="non-scaling-stroke"`);
  addLine(out, 130, trackY.inv, 1430, trackY.inv, `stroke="#d0d0d0" stroke-width="1" vector-effect="non-scaling-stroke"`);

  for (const y of Object.values(trackY)) {
    addCollapsedBlock(out, leftX, y, 85, "WIPI2-CCZ1B");
    addCollapsedBlock(out, rightFlankX, y, 120, "CCZ1B-NXPH1");
  }

  addPms2Gene(out, g, eventConf, { name: "RSPH10B", rev: true }, blockX - 74, trackY.ref, 55, true);
  addPms2Gene(out, g, eventConf, { name: "RSPH10B", rev: true }, blockX - 74, trackY.inv, 55, true);
  const refPos = drawPms2GeneSeries(out, g, eventConf, refEntries, blockX, trackY.ref, geneWidth, gap, highlightGenes);
  const invPos = drawPms2GeneSeries(out, g, eventConf, invEntries, blockX, trackY.inv, geneWidth, gap, highlightGenes);
  addPms2Gene(out, g, eventConf, { name: "RSPH10B2", rev: false }, blockEnd + 16, trackY.ref, 58, true);
  addPms2Gene(out, g, eventConf, { name: "RSPH10B2", rev: false }, blockEnd + 16, trackY.inv, 58, true);

  addBracket(out, blockX, blockEnd, 70, "17-gene block", "#555555");
  addBracket(out, blockX, blockEnd, 190, "same genes inverted", "#0072B2");
  addLine(
    out,
    blockX - 5,
    178,
    blockEnd + 5,
    178,
    `stroke="#0072B2" stroke-width="2" stroke-dasharray="6 4" vector-effect="non-scaling-stroke"`
  );
  addText(
    out,
    (blockX + blockEnd) / 2,
    169,
    "large inversion mediated by PMS2/pseudogene segmental duplication",
    `text-anchor="middle" font-family="monospace" font-size="10pt" fill="#0072B2"`
  );

  addBracket(out, blockX - 74, blockX + geneWidth, 295, "PMS2-side duplicated copy", "#6A51A3");
  addBracket(out, blockEnd - geneWidth + 8, blockEnd + 74, 295, "pseudogene-side duplicated copy", "#6A51A3");

  const refPms2 = refPos.get("PMS2");
  const invPms2 = invPos.get("PMS2");
  const leftCopyX = blockX - 45;
  const rightCopyX = blockEnd + 45;
  addPath(
    out,
    `M ${refPms2.x + refPms2.width / 2} ${trackY.ref + 30} C ${leftCopyX + 120} 315 ${rightCopyX - 120} 315 ${rightCopyX} ${trackY.inv + 30}`,
    `fill="none" stroke="#D55E00" stroke-width="1.5" stroke-dasharray="5 4" vector-effect="non-scaling-stroke"`
  );
  addPath(
    out,
    `M ${rightCopyX} ${trackY.ref + 35} C ${rightCopyX - 120} 340 ${leftCopyX + 120} 340 ${invPms2.x + invPms2.width / 2} ${trackY.inv + 35}`,
    `fill="none" stroke="#D55E00" stroke-width="1.5" stroke-dasharray="5 4" vector-effect="non-scaling-stroke"`
  );
  addText(
    out,
    (blockX + blockEnd) / 2,
    336,
    "minimap2 can prefer contiguous alignments, swapping PMS2 and pseudogene mappings",
    `text-anchor="middle" font-family="monospace" font-size="10pt" fill="#D55E00"`
  );

  out.push(`</svg>`);
  return out.join("\n") + "\n";
}

function renderGraphSvg(conf, g, options) {
  const [sub, pos] = calPos(conf, g);
  for (let i = 0; i < pos.length; ++i) {
    pos[i].cx_st = conf.xskip + pos[i].start;
    pos[i].cx_en = pos[i].cx_st + pos[i].len;
    pos[i].cy = conf.yskip * (pos[i].level + 1);
  }

  let maxW = 0;
  let maxH = 0;
  for (let i = 0; i < pos.length; ++i) {
    maxW = Math.max(maxW, pos[i].cx_en + conf.xskip);
    maxH = Math.max(maxH, pos[i].cy + conf.yskip);
  }

  const r2c = rankToColor(g);
  const out = [];
  out.push(`<?xml version="1.0" encoding="UTF-8"?>`);
  out.push(
    `<svg xmlns="http://www.w3.org/2000/svg" width="${maxW}" height="${maxH}" ` +
      `viewBox="0 0 ${maxW} ${maxH}">`
  );
  if (options.title) out.push(`  <title>${xmlEscape(options.title)}</title>`);
  if (options.background !== "none") {
    out.push(`  <rect width="100%" height="100%" fill="${xmlEscape(options.background)}" />`);
  }

  out.push(`  <g id="edges" opacity="0.8" fill="none">`);
  for (let i = 0; i < pos.length; ++i) {
    for (let j = 0; j < sub.v[i].n; ++j) {
      const k = sub.a[sub.v[i].off + j].i;
      const r = sub.a[sub.v[i].off + j].rank;
      const color = r >= 0 && r2c[r] != null ? r2c[r] : "#A0A0A0";
      out.push(
        `    <line x1="${pos[i].cx_en}" y1="${pos[i].cy}" ` +
          `x2="${pos[k].cx_st}" y2="${pos[k].cy}" stroke="${color}" />`
      );
    }
  }
  out.push(`  </g>`);

  out.push(`  <g id="genes">`);
  for (let i = 0; i < pos.length; ++i) {
    const s = g.seg[sub.v[i].v >> 1];
    const label = conf.label === "name" ? s.name : s.len;
    const stroke = s.rank >= 0 && r2c[s.rank] != null ? r2c[s.rank] : null;
    svgArrow(
      out,
      pos[i].cx_st,
      pos[i].cy,
      pos[i].len,
      conf.h_arrow,
      Boolean(sub.v[i].v & 1),
      label,
      conf.font_size,
      stroke,
      s.color
    );
  }
  out.push(`  </g>`);
  out.push(`</svg>`);
  return out.join("\n") + "\n";
}

function main() {
  const args = parseArgs(process.argv);
  const conf = plotConf(args);
  const gfa = fs.readFileSync(args.input, "utf8");
  const g = parseGfa(gfa);
  const options = {
    background: args.background,
    title: args.title,
    referenceWalk: args.referenceWalk,
    highlightWalks: args.highlightWalks,
    fadeOpacity: args.fadeOpacity,
  };
  const svg =
    args.layout === "reference-walk"
      ? renderReferenceWalkSvg(conf, g, options)
      : args.layout === "event-tracks"
        ? renderEventTracksSvg(conf, g, options)
        : args.layout === "pms2-events"
          ? renderPms2EventsSvg(conf, g, options)
      : renderGraphSvg(conf, g, options);

  fs.mkdirSync(path.dirname(args.output), { recursive: true });
  fs.writeFileSync(args.output, svg);
}

try {
  main();
} catch (err) {
  console.error(`pangene_gfa_to_svg.js: ${err.message}`);
  process.exit(1);
}
