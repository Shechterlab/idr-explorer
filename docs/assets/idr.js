/* IDR Explorer (Shechter Lab). Static, client-side app over docs/data/, which Analysis/scripts/build_idr_explorer.py
   (Shechterlab/MicroDomains) writes from idr_universe.py and idr_overview_stats.py.

   The Overview is the Methylarginine Site Explorer's paradigm with that site's own code (ptm-core.js, generated from
   its explorer.js by tools/sync-ptm-core.mjs): a selection narrows the proteins and the evidence, and every tile's test
   runs again on it (idr-stats.js), drawn with its scorecard tiles (aHome, aForest), charts (hbars) and exports. The
   per-protein files use its 64 shards and hash, and its domain rows. */
(() => {
  'use strict';
  const C = window.PTMCore;
  const $ = (s, el = document) => el.querySelector(s);
  const $$ = (s, el = document) => Array.from(el.querySelectorAll(s));
  const esc = C.esc, fmt = C.fmt, fmtP = C.fmtP, fx = C.fx;
  const PADJ = C.PADJ;
  const SHARDS = 64;
  const shardOf = (acc) => { let h = 0; for (const c of acc) h = (h * 31 + c.charCodeAt(0)) % 1000003; return h % SHARDS; };
  const PTM_SITE = 'https://shechterlab.github.io/PTM_Pipeline/';
  const GH = { site: 'https://github.com/Shechterlab/idr-explorer/blob/main/', ptm: 'https://github.com/Shechterlab/PTM_Pipeline/blob/main/' };
  const PAGE = 100;
  const ASSET_V = (() => { try { return new URL(document.currentScript.src).searchParams.get('v') || ''; } catch (e) { return ''; } })();

  // ---------------------------------------------------------------- data
  const cache = new Map();
  const getJSON = (path) => {
    if (!cache.has(path)) cache.set(path, fetch(path + (ASSET_V ? '?v=' + ASSET_V : '')).then((r) => { if (!r.ok) throw new Error(path + ': HTTP ' + r.status); return r.json(); }));
    return cache.get(path);
  };
  let B = null, OV = null, PR = null, FL = null, ID = null, SI = null, RS = null, ST = null;
  const byAcc = new Map(), byGene = new Map();
  const loadBase = async () => {
    if (!B) [B, OV] = await Promise.all([getJSON('data/build.json'), getJSON('data/overview.json')]);
  };
  const loadProteins = async () => {
    if (!PR) {
      PR = await getJSON('data/proteins.json');
      PR.acc.forEach((a, i) => byAcc.set(a, i));
      PR.gene.forEach((g, i) => { const k = (g || '').toUpperCase(); if (k && !byGene.has(k)) byGene.set(k, i); });
    }
    return PR;
  };
  const loadFlanks = async () => (FL = FL || await getJSON('data/flanks.json'));
  const loadIdrs = async () => (ID = ID || await getJSON('data/idrs.json'));
  const loadShard = (acc) => getJSON('data/p/' + String(shardOf(acc)).padStart(2, '0') + '.json');
  // the data the Overview's tests need for a selection (about 12 MB): loaded on the first selection
  const loadStats = async () => {
    if (!ST) {
      [SI, RS] = await Promise.all([getJSON('data/sites.json'), getJSON('data/residues.json'), loadProteins(), loadIdrs()]);
      ST = window.IDRStats.make(C, { overview: OV, sites: SI, residues: RS, proteins: PR, idrs: ID });
    }
    return ST;
  };

  // ---------------------------------------------------------------- vocabulary (overview.json)
  const ctxLabel = (c) => (OV.contexts.find((x) => x[0] === c) || [c, c])[1];
  const ctxDef = (c) => (OV.contexts.find((x) => x[0] === c) || [c, c, ''])[2];
  const kindLabel = (k) => (OV.kinds.find((x) => x[0] === k) || [k, k])[1];
  const CTX_ORDER = () => OV.contexts.map((c) => c[0]);
  const srcBit = (k) => (OV.sources.find((s) => s[0] === k) || [k, 0])[1];
  const srcNames = (bits) => OV.sources.filter((s) => bits & s[1]).map((s) => s[2]).join(', ');
  const CONF = ['', 'Bronze', 'Silver', 'Gold'];
  const tableRows = (t) => t.rows.map((r) => Object.fromEntries(t.cols.map((c, i) => [c, r[i]])));
  const l2 = (v) => (v != null && v > 0 && Number.isFinite(v) ? +Math.log2(v).toFixed(3) : null);
  const CHEM_SHORT = ['Unmarked', 'R-type', 'K-type', 'Hybrid'];
  const FLAGS = [[1, 'TF'], [2, 'RBD'], [8, 'chaperone'], [16, 'HMG'], [4, 'histone-like']];
  const flagText = (f) => FLAGS.filter(([b]) => f & b).map(([, t]) => t).join(', ');
  const geneOf = (p) => PR.gene[p] || PR.acc[p];
  const protLink = (p, extra) => '<a href="#' + esc(PR.acc[p]) + (extra ? '?' + extra : '') + '">' + esc(geneOf(p)) + '</a>';
  const z = (v) => (v == null ? null : v / 100);
  const fz = (v) => (v == null ? '' : (v >= 0 ? '+' : '−') + Math.abs(v / 100).toFixed(2));
  const zChip = (v, pole) => {
    if (v == null) return '<span class="muted">–</span>';
    const a = Math.max(0, Math.min(1, (v / 100) / 3));
    return '<span class="z" style="background:color-mix(in srgb, var(--' + pole + ') ' + Math.round(a * 55) + '%, transparent)">' + fz(v) + '</span>';
  };

  // ---------------------------------------------------------------- selection (the explorer's F and R)
  // Proteins: gene list, regions they have (all must hold), classes (any of those ticked). Evidence: the sources a site
  // must come from (any of them), and the lowest PTMeXchange grade for a site PTMeXchange alone reports.
  const HAS_LABEL = { acidic_pole: 'an acidic-pole window', basic_pole: 'a basic-pole window', blocky: 'a blocky IDR', acidic_stretch: 'an acidic stretch',
    basic_patch: 'a basic patch', rg_cluster: 'an RG/RGG cluster' };
  const F = { genes: '', has: [], cls: [], src: null, conf: 0 };
  const R = { prots: [], sites: [], scoped: false, sel: null };
  let SELKEY = '';
  const readSel = (params) => {
    F.genes = params.get('genes') || '';
    F.has = (params.get('has') || '').split(',').filter(Boolean);
    F.cls = (params.get('cls') || '').split(',').filter(Boolean);
    F.src = params.has('src') ? +params.get('src') : null;
    F.conf = +(params.get('conf') || 0);
  };
  const selParams = () => {
    const p = new URLSearchParams();
    if (F.genes) p.set('genes', F.genes);
    if (F.has.length) p.set('has', F.has.join(','));
    if (F.cls.length) p.set('cls', F.cls.join(','));
    if (F.src != null) p.set('src', F.src);
    if (F.conf) p.set('conf', F.conf);
    return p;
  };
  const allSrc = () => OV.sources.reduce((b, s) => b | s[1], 0);
  const geneTokens = () => F.genes.split(/[\s,;]+/).map((t) => t.trim().toUpperCase()).filter(Boolean);
  function constraints() {
    if (!OV) return [];
    const out = [];
    if (F.genes) out.push('proteins: ' + geneTokens().join(', '));
    for (const h of F.has) out.push('with ' + (HAS_LABEL[h] || h));
    if (F.cls.length) out.push(F.cls.map((k) => (OV.classes.find((c) => c[0] === k) || [k, k])[1]).join(' or '));
    if (F.src != null && F.src !== allSrc()) out.push('sites from ' + srcNames(F.src));
    if (F.conf) out.push('PTMeXchange-only sites ' + CONF[F.conf] + ' or better');
    return out;
  }
  const selectionText = () => constraints().join('; ') || 'all proteins and sites';
  // the selection's proteins (indices, or null for all) and the evidence test of a site
  function protSelection() {
    if (!F.genes && !F.has.length && !F.cls.length) return null;
    const toks = geneTokens(), want = new Set();
    for (const t of toks) {
      if (byAcc.has(t)) want.add(byAcc.get(t));
      else if (byGene.has(t)) want.add(byGene.get(t));
      else if (toks.length === 1) PR.gene.forEach((g, i) => { if ((g || '').toUpperCase().startsWith(t)) want.add(i); });
    }
    const hasBits = F.has.map((h) => 1 << OV.has.indexOf(h)), clsBits = F.cls.map((k) => 1 << OV.classes.findIndex((c) => c[0] === k));
    const out = new Set();
    for (let i = 0; i < PR.acc.length; i++) {
      if (toks.length && !want.has(i)) continue;
      if (hasBits.some((b) => !(PR.has[i] & b))) continue;
      if (clsBits.length && !clsBits.some((b) => PR.cls[i] & b)) continue;
      out.add(i);
    }
    return out;
  }
  function evidenceKeep() {
    const src = F.src == null ? allSrc() : F.src, ptmx = srcBit('ptmx');
    if (src === allSrc() && !F.conf) return null;
    return (i) => (SI.src[i] & src) !== 0 && !(F.conf && SI.src[i] === ptmx && SI.conf[i] < F.conf);
  }
  async function computeSelection() {
    const key = selParams().toString();
    R.scoped = constraints().length > 0;
    if (!R.scoped) { R.sel = null; SELKEY = key; R.prots = []; R.sites = []; return; }
    if (key === SELKEY && R.sel) return;
    await loadProteins();
    await loadStats();
    const prots = protSelection(), keep = evidenceKeep();
    R.sel = { prots, keep };
    R.prots = prots ? [...prots] : PR.acc.map((_, i) => i);
    const sites = [];
    for (let i = 0; i < SI.p.length; i++) if ((!prots || prots.has(SI.p[i])) && (!keep || keep(i))) sites.push(i);
    R.sites = sites;
    SELKEY = key;
  }
  const store = { get(k) { try { return localStorage.getItem('idr.' + k); } catch (e) { return null; } }, set(k, v) { try { if (v == null) localStorage.removeItem('idr.' + k); else localStorage.setItem('idr.' + k, v); } catch (e) { /* storage unavailable */ } } };
  const docLink = (key) => ' <a class="linkish doc-link" href="#definitions/' + esc(key) + '">Method</a>';
  const UI = { charts: {}, view: 'overview' };
  const sizeOf = (n) => ({ get length() { return n(); } });
  C.use({ M: {}, UI, R, F, S: sizeOf(() => (B ? B.counts.sites : 0)), P: sizeOf(() => (B ? B.counts.proteins : 0)), constraints, selectionText, docLink,
    creditLines: () => ['IDR Explorer · shechterlab.github.io/idr-explorer', 'Shechter Lab · cite the data sources (Sources page)'], evHoverCard: () => null, store });

  // ---------------------------------------------------------------- routing
  const parseHash = () => {
    let h = location.hash.slice(1);
    try { h = decodeURIComponent(h); } catch (e) { /* keep */ }
    const [path, qs] = h.split('?');
    return { path: path || 'overview', params: new URLSearchParams(qs || '') };
  };
  const setParams = (path, params) => { const s = params.toString(); history.replaceState(null, '', '#' + path + (s ? '?' + s : '')); };
  const view = () => $('#view');
  const setTab = (name) => $$('.tabs a').forEach((a) => {
    if (a.dataset.view === name) a.setAttribute('aria-current', 'page'); else a.removeAttribute('aria-current');
    const s = selParams().toString();
    a.setAttribute('href', '#' + a.dataset.view + (s ? '?' + s : ''));
  });
  let routeToken = 0;
  async function route() {
    const token = ++routeToken;
    const { path, params } = parseHash();
    const [base, sub] = path.split('/');
    const views = { overview: overviewView, flanks: flanksView, idrs: idrsView, definitions: definitionsView, about: aboutView };
    const v = view();
    v.innerHTML = '<p class="loading">Loading…</p>';
    try {
      await loadBase();
      if (token !== routeToken) return;
      readSel(params);
      setTab(views[base] ? base : '');
      document.title = views[base] ? 'IDR Explorer · Shechter Lab' : document.title;
      if (views[base]) await views[base](params, token, sub);
      else await proteinView(path, params, token);
      if (sub && base === 'definitions') { const t = document.getElementById('def-' + sub); if (t) t.scrollIntoView(); }
    } catch (e) {
      if (token === routeToken) v.innerHTML = '<div class="panel"><b>Could not load this view.</b> <span class="muted">' + esc(e.message) + '</span></div>';
      console.error(e);
    }
  }
  window.addEventListener('hashchange', route);

  // ---------------------------------------------------------------- tooltips and exports (the explorer's)
  document.addEventListener('mouseover', (e) => { const t = e.target.closest && e.target.closest('[data-tip]'); if (t) C.showTip(e, t.dataset.tip); });
  document.addEventListener('mousemove', (e) => { if (e.target.closest && e.target.closest('[data-tip]')) C.moveTip(e); });
  document.addEventListener('mouseout', (e) => { const t = e.target.closest && e.target.closest('[data-tip]'); if (t && !t.contains(e.relatedTarget)) C.hideTip(); });
  function bindExports(root) {
    $$('[data-exp][data-chart]', root).forEach((b) => b.addEventListener('click', () => {
      const spec = UI.charts[b.dataset.chart];
      const name = 'idr_' + spec.id.replace(/^ch-/, '');
      if (b.dataset.exp === 'tsv') { C.download(name + '.tsv', spec.figTsv ? spec.figTsv() : C.chartTsv(spec)); return; }
      C.exportFigure(b.dataset.exp, spec.figSvg ? C.figExportSvg(spec) : C.chartSvg(spec), name);
    }));
  }
  // a chart card as the explorer draws one (chartCard without its clusterProfiler note): title, subtitle, exports, bars
  function card(spec, extra) {
    UI.charts[spec.id] = spec;
    return '<div class="card chart-card" id="' + spec.id + '"><div class="card-head"><div><h3>' + esc(spec.title) + '</h3><p class="sub">' + spec.subHtml + '</p></div>' + C.expButtons(spec.id) + '</div>' +
      (spec.tools ? '<div class="chart-tools">' + spec.tools + '</div>' : '') + (spec.figSvg ? spec.figSvg() : C.hbars(spec)) + (spec.footHtml || '') + (extra || '') + '</div>';
  }

  // ---------------------------------------------------------------- shared table
  function table(host, opts) {
    const { caption, columns, rows, sortKey, desc, rowHeader = 0, csvName } = opts;
    let { page } = opts;
    const onState = opts.onState || ((o) => table(host, { ...opts, sortKey: o.sort ?? sortKey, desc: o.desc ?? desc, page: o.page ?? page }));
    const order = rows.slice();
    const col = columns.find((c) => c.key === sortKey) || columns[0];
    const sv = col.sort || col.get;
    order.sort((a, b) => {
      const x = sv(a), y = sv(b);
      if (x == null && y == null) return 0;
      if (x == null) return 1;
      if (y == null) return -1;
      return (x < y ? -1 : x > y ? 1 : 0) * (desc ? -1 : 1);
    });
    const pages = Math.max(1, Math.ceil(order.length / PAGE));
    page = Math.min(Math.max(1, page || 1), pages);
    const slice = order.slice((page - 1) * PAGE, page * PAGE);
    const id = 'cap-' + Math.random().toString(36).slice(2, 8);
    const head = columns.map((c) => {
      const sorted = c.key === col.key ? (desc ? 'descending' : 'ascending') : null;
      return '<th scope="col"' + (sorted ? ' aria-sort="' + sorted + '"' : '') + (c.tip ? ' data-tip="' + esc(c.tip) + '"' : '') + '><button type="button" data-sort="' + c.key + '">' + esc(c.label) + '</button></th>';
    }).join('');
    const body = slice.map((i) => '<tr>' + columns.map((c, k) => {
      const cell = c.html ? c.html(i) : esc(c.get(i) ?? '');
      return k === rowHeader ? '<th scope="row" data-label="' + esc(c.label) + '">' + cell + '</th>' : '<td' + (c.num ? ' class="num"' : '') + ' data-label="' + esc(c.label) + '">' + cell + '</td>';
    }).join('') + '</tr>').join('');
    host.innerHTML = '<div class="table-wrapper" tabindex="0" role="region" aria-labelledby="' + id + '"><table class="data"><caption id="' + id + '">' + esc(caption) +
      '</caption><thead><tr>' + head + '</tr></thead><tbody>' + body + '</tbody></table></div>' +
      '<div class="pager"><button class="btn" type="button" data-page="' + (page - 1) + '"' + (page <= 1 ? ' disabled' : '') + '>Previous</button>' +
      '<span>Page ' + page + ' of ' + fmt(pages) + ' · ' + fmt(order.length) + ' rows</span>' +
      '<button class="btn" type="button" data-page="' + (page + 1) + '"' + (page >= pages ? ' disabled' : '') + '>Next</button>' +
      (csvName ? '<button class="btn" type="button" data-csv="1">Download these ' + fmt(order.length) + ' rows (TSV)</button>' : '') + '</div>';
    $$('[data-sort]', host).forEach((b) => b.addEventListener('click', () => { const k = b.dataset.sort; onState({ sort: k, desc: k === col.key ? !desc : true, page: 1 }); }));
    $$('[data-page]', host).forEach((b) => b.addEventListener('click', () => onState({ page: +b.dataset.page })));
    const csv = $('[data-csv]', host);
    if (csv) csv.addEventListener('click', () => C.download(csvName, [columns.map((c) => c.label).join('\t')].concat(order.map((i) => columns.map((c) => String(c.get(i) ?? '')).join('\t'))).join('\n') + '\n'));
  }

  // A binned density of two columns (all rows in grey, the selection in the accent), as SVG.
  function density({ xs, ys, sel, xLabel, yLabel, xr, yr, bins = 56 }) {
    const W = 420, H = 340, L = 46, Rm = 10, T = 10, Bm = 40, pw = W - L - Rm, ph = H - T - Bm;
    const bx = (x) => Math.min(bins - 1, Math.max(0, Math.floor((x - xr[0]) / (xr[1] - xr[0]) * bins)));
    const by = (y) => Math.min(bins - 1, Math.max(0, Math.floor((y - yr[0]) / (yr[1] - yr[0]) * bins)));
    const grid = (idx) => { const g = new Uint32Array(bins * bins); for (const i of idx) { const x = xs[i], y = ys[i]; if (x == null || y == null) continue; g[by(y) * bins + bx(x)]++; } return g; };
    const rect = (g, color, scale, what) => {
      let max = 0; for (const v of g) if (v > max) max = v;
      let out = '';
      const cw = pw / bins, ch = ph / bins;
      for (let k = 0; k < g.length; k++) {
        if (!g[k]) continue;
        const a = 0.12 + 0.88 * Math.log1p(g[k]) / Math.log1p(max || 1), gx = k % bins, gy = Math.floor(k / bins);
        out += '<rect x="' + (L + gx * cw).toFixed(1) + '" y="' + (T + ph - (gy + 1) * ch).toFixed(1) + '" width="' + (cw + 0.3).toFixed(1) + '" height="' + (ch + 0.3).toFixed(1) +
          '" style="fill:' + color + ';fill-opacity:' + (a * scale).toFixed(2) + '" data-tip="' + esc(fmt(g[k]) + ' ' + what) + '"/>';
      }
      return out;
    };
    const sx = (x) => L + (x - xr[0]) / (xr[1] - xr[0]) * pw, sy = (y) => T + ph - (y - yr[0]) / (yr[1] - yr[0]) * ph;
    const ticks = (r) => { const out = []; for (let t = Math.ceil(r[0]); t <= Math.floor(r[1]); t++) out.push(t); return out; };
    let svg = '<svg class="plane" viewBox="0 0 ' + W + ' ' + H + '" role="img" aria-label="' + esc(yLabel + ' against ' + xLabel) + '"><rect x="' + L + '" y="' + T + '" width="' + pw + '" height="' + ph + '" style="fill:var(--surface-2)"/>';
    svg += rect(grid(xs.keys()), 'var(--muted)', sel ? 0.55 : 1, 'rows');
    if (sel) svg += rect(grid(sel), 'var(--accent)', 1, 'selected rows');
    for (const t of ticks(xr)) svg += '<line x1="' + sx(t) + '" x2="' + sx(t) + '" y1="' + (T + ph) + '" y2="' + (T + ph + 4) + '" style="stroke:var(--muted)"/><text x="' + sx(t) + '" y="' + (T + ph + 16) + '" text-anchor="middle" font-size="12">' + t + '</text>';
    for (const t of ticks(yr)) svg += '<line x1="' + (L - 4) + '" x2="' + L + '" y1="' + sy(t) + '" y2="' + sy(t) + '" style="stroke:var(--muted)"/><text x="' + (L - 7) + '" y="' + (sy(t) + 4) + '" text-anchor="end" font-size="12">' + t + '</text>';
    if (xr[0] < 0 && xr[1] > 0) svg += '<line x1="' + sx(0) + '" x2="' + sx(0) + '" y1="' + T + '" y2="' + (T + ph) + '" style="stroke:var(--line-strong);stroke-dasharray:3 3"/>';
    if (yr[0] < 0 && yr[1] > 0) svg += '<line x1="' + L + '" x2="' + (L + pw) + '" y1="' + sy(0) + '" y2="' + sy(0) + '" style="stroke:var(--line-strong);stroke-dasharray:3 3"/>';
    svg += '<text x="' + (L + pw / 2) + '" y="' + (H - 6) + '" text-anchor="middle" font-size="12">' + esc(xLabel) + '</text><text transform="translate(14 ' + (T + ph / 2) + ') rotate(-90)" text-anchor="middle" font-size="12">' + esc(yLabel) + '</text>';
    return svg + '</svg>';
  }

  // Horizontal box plots from quantiles [5th, 25th, median, 75th, 95th percentile, n] (the explorer's dmBoxSvg layout):
  // box the middle half, whiskers the 5th to 95th percentile, bar the median. rows: [label, q, cls, tip].
  function boxSvg(rows, axis, f) {
    const W = 520, ML = 190, MR = 14, RH = 24, MT = 6, MB = 34, H = MT + rows.length * RH + MB;
    const vals = rows.flatMap((r) => (r[1] ? [r[1][0], r[1][4]] : []));
    let lo = Math.min(0, ...vals), hi = Math.max(...vals);
    const pad = (hi - lo) * 0.04 || 0.05;
    hi += pad;
    const x = (v) => ML + ((v - lo) / (hi - lo)) * (W - ML - MR);
    const step = [0.01, 0.02, 0.05, 0.1, 0.2, 0.25, 0.5, 1, 2, 5, 10, 20, 25, 50].find((st) => (hi - lo) / st <= 6) || 100;
    let s = '';
    rows.forEach(([label, q, cls, tp], k) => {
      const cy = MT + k * RH + RH / 2;
      s += '<text class="bx-lab" x="' + (ML - 8) + '" y="' + (cy + 4).toFixed(1) + '" text-anchor="end">' + esc(label) + (q ? ' (' + fmt(q[5]) + ')' : '') + '</text>';
      if (!q) { s += '<text class="bx-tick" x="' + ML + '" y="' + (cy + 4) + '">none</text>'; return; }
      const t = ' data-tip="' + esc((tp ? tp + '<br>' : '') + 'median ' + f(q[2]) + '; middle half ' + f(q[1]) + ' to ' + f(q[3]) + '; 5th–95th percentile ' + f(q[0]) + ' to ' + f(q[4]) + '; ' + fmt(q[5]) + ' values') + '"';
      s += '<g class="bx ' + (cls || '') + '"' + t + '><line x1="' + x(q[0]).toFixed(1) + '" x2="' + x(q[4]).toFixed(1) + '" y1="' + cy + '" y2="' + cy + '"/>' +
        '<rect x="' + x(q[1]).toFixed(1) + '" y="' + (cy - 7) + '" width="' + Math.max(1.5, x(q[3]) - x(q[1])).toFixed(1) + '" height="14" rx="2"/>' +
        '<line class="bx-med" x1="' + x(q[2]).toFixed(1) + '" x2="' + x(q[2]).toFixed(1) + '" y1="' + (cy - 8) + '" y2="' + (cy + 8) + '"/>' +
        '<rect x="' + x(q[0]).toFixed(1) + '" y="' + (cy - 10) + '" width="' + Math.max(4, x(q[4]) - x(q[0])).toFixed(1) + '" height="20" fill="transparent"/></g>';
    });
    s += '<line class="bx-axis" x1="' + ML + '" x2="' + (W - MR) + '" y1="' + (H - MB) + '" y2="' + (H - MB) + '"/>';
    for (let v = Math.ceil(lo / step) * step; v <= hi + 1e-12; v += step) s += '<line class="bx-axis" x1="' + x(v).toFixed(1) + '" x2="' + x(v).toFixed(1) + '" y1="' + (H - MB) + '" y2="' + (H - MB + 4) + '"/><text class="bx-tick" x="' + x(v).toFixed(1) + '" y="' + (H - MB + 15) + '" text-anchor="middle">' + f(v) + '</text>';
    s += '<text class="bx-axt" x="' + ML + '" y="' + (H - 3) + '">' + esc(axis) + '</text>';
    return '<svg class="bx-plot" viewBox="0 0 ' + W + ' ' + H + '" role="img" aria-label="' + esc(axis) + '">' + s + '</svg>';
  }
  const BOX_CSS = (pal) => '.bx line{stroke:' + pal.ink + ';stroke-width:1.2}.bx rect{fill:#b9dde6;stroke:#0f3b47}.bx.in rect{fill:#f8dede;stroke:#d64545}.bx-med{stroke:' + pal.ink + ';stroke-width:2}' +
    '.bx-lab,.bx-tick,.bx-axt{font-size:11px;fill:' + pal.muted + '}.bx-lab{fill:' + pal.ink + '}.bx-axis{stroke:' + pal.muted + '}';

  // ---------------------------------------------------------------- the Overview's tests as charts
  // A row of a chart (hbars) and a tile (aForest): value log2 of the effect, q, and the counts behind it.
  const SIG = (q) => q != null && q <= PADJ;
  function ctxRows(ctx, sel) {
    const kinds = ctx === 'basophilic' ? OV.kinds.filter((k) => 'RK'.includes(k[2])) : OV.kinds;
    const whole = new Map(tableRows(OV.residue_tests).filter((r) => r.context === ctx).map((r) => [r.kind, r]));
    const val = (x) => (x.tested ? l2(x.mh_or) : null);
    let rows = kinds.map(([k, label, res]) => {
      const w = whole.get(k) || {};
      return { key: k, label, note: res, value: val(w), q: w.q_in_context, sig: SIG(w.q_in_context), w, nums: [w.sites_in, w.sites], extra: fmt(w.sites_in) + '/' + fmt(w.sites) };
    }).filter((r) => r.w.sites);
    if (sel) {
      const res = rows.map((r) => ST.residueTest(r.key, ctx, sel));
      C.bhAdjust(res.filter((x) => x.tested), (x) => x.mh_p);
      rows.forEach((r, i) => {
        const x = res[i];
        r.ref = r.value; r.refQ = r.q; r.s = x;
        r.value = val(x); r.q = x.tested ? x.q : null; r.sig = SIG(r.q); r.nums = [x.sites_in, x.sites]; r.extra = fmt(x.sites_in) + '/' + fmt(x.sites);
      });
    }
    for (const r of rows) {
      const x = r.s || r.w;
      const few = x.tested ? '' : '<br><i>Too few to test: ' + fmt(x.sites_in + x.other_in) + ' residues in this context and ' + fmt(x.mh_strata) + ' informative proteins (a test needs ' + OV.settings.min_in_context + ' and ' + OV.settings.min_strata + ')</i>';
      r.title = '<b>' + esc(r.label) + ' in ' + esc(ctxLabel(ctx).toLowerCase()) + '</b>' + few + '<br>' + fmt(x.sites_in) + ' of ' + fmt(x.sites) + ' sites, against ' + fmt(x.other_in) + ' of ' + fmt(x.other) +
        ' unmodified ' + esc(r.note.split('').join('/')) + ' in the same ' + fmt(x.proteins) + ' proteins<br>within protein (Mantel–Haenszel): OR ' + C.fmtOr(x.mh_or, x.mh_lo, x.mh_hi) +
        ', q ' + fmtP(r.q) + ' · ' + fmt(x.mh_strata) + ' informative proteins<br>pooled: OR ' + (x.pooled_or == null ? '—' : x.pooled_or.toFixed(2)) + ', Fisher p ' + fmtP(x.pooled_p) +
        (r.s ? '<br><span class="muted">Whole dataset: OR ' + C.fmtOr(r.w.mh_or, r.w.mh_lo, r.w.mh_hi) + ', q ' + fmtP(r.refQ) + ' (tick)</span>' : '');
      r.main = { or: x.mh_or, lo: x.mh_lo, hi: x.mh_hi, q: r.q, sig: r.sig, used: x.mh_strata };
    }
    return rows;
  }
  function selTestOf(sel, unit, text, method) {
    return { n: R.sites.length, N: B.counts.sites, unit, what: selectionText(), text, method };
  }
  function ctxSpec(ctx) {
    const rows = ctxRows(ctx, R.sel);
    rows.sort((a, b) => (b.value ?? -99) - (a.value ?? -99));
    const spec = { id: 'ch-' + ctx, title: ctxLabel(ctx), rows, extraCol: true, heads: ['Modification', 'in / sites', '', 'log2 OR'], extraName: 'sites_in_context_of_sites',
      subHtml: esc(ctxDef(ctx)) + '. Each modification’s sites against the unmodified residues of the same type in the same proteins: odds ratio within protein (Mantel–Haenszel), BH across the modifications.' + docLink('tests'),
      subText: ctxLabel(ctx) + ': ' + ctxDef(ctx) + '. Sites of each modification against the unmodified residues of the same type in the same proteins, within protein (Mantel–Haenszel), BH across modifications.',
      axis: 'log2 odds ratio, within protein (sites against unmodified residues of the same type)', valueLabel: 'log2_mh_odds_ratio',
      tsvCols: [['proteins', (r) => (r.s || r.w).proteins], ['mh_strata', (r) => (r.s || r.w).mh_strata], ['mh_or', (r) => (r.s || r.w).mh_or], ['mh_ci_low', (r) => (r.s || r.w).mh_lo],
        ['mh_ci_high', (r) => (r.s || r.w).mh_hi], ['mh_p', (r) => (r.s || r.w).mh_p], ['pooled_or', (r) => (r.s || r.w).pooled_or], ['pooled_fisher_p', (r) => (r.s || r.w).pooled_p]] };
    if (R.sel) spec.selTest = selTestOf(R.sel, 'sites', fmt(R.sites.length) + ' sites of your ' + fmt(R.prots.length) + ' proteins, against the unmodified residues of the same type in the same proteins',
      'Mantel–Haenszel odds ratio with each protein a stratum, BH across the modifications');
    return spec;
  }
  function kinaseSpec(ctx) {
    const W = tableRows(OV.kinase_tests).filter((r) => r.context === ctx);
    const rows = OV.kinase_rows.map(([key, label, cls], b) => {
      const w = W.find((x) => x.row === key) || {};
      return { key, label, note: cls, value: w.tested ? l2(w.mh_or) : null, q: w.q, sig: SIG(w.q), w, extra: fmt(w.in_context_match) + '/' + fmt(w.in_context), b };
    });
    if (R.sel) {
      const res = rows.map((r) => ST.kinaseTest(ctx, r.b, R.sel));
      C.bhAdjust(res.filter((x) => x.tested), (x) => x.mh_p);
      rows.forEach((r, i) => { const x = res[i]; r.ref = r.value; r.refQ = r.q; r.s = x; r.value = x.tested ? l2(x.mh_or) : null; r.q = x.tested ? x.q : null; r.sig = SIG(r.q); r.extra = fmt(x.in_context_match) + '/' + fmt(x.in_context); });
    }
    for (const r of rows) {
      const x = r.s || r.w;
      r.title = '<b>' + esc(r.label) + '</b> (' + esc(r.note) + ') motif matches (percentile ≥ ' + OV.kinase_match + ')<br>' + fmt(x.in_context_match) + ' of ' + fmt(x.in_context) + ' phosphosites in ' + esc(ctxLabel(ctx).toLowerCase()) +
        ', against ' + fmt(x.outside_match) + ' of ' + fmt(x.outside) + ' elsewhere in the same proteins<br>within protein: OR ' + C.fmtOr(x.mh_or, x.mh_lo, x.mh_hi) + ', q ' + fmtP(r.q) +
        (r.s ? '<br><span class="muted">Whole dataset: OR ' + C.fmtOr(r.w.mh_or, r.w.mh_lo, r.w.mh_hi) + ', q ' + fmtP(r.refQ) + ' (tick)</span>' : '');
      r.main = { or: x.mh_or, lo: x.mh_lo, hi: x.mh_hi, q: r.q, sig: r.sig, used: x.mh_strata };
    }
    const spec = { id: 'ch-kin-' + ctx, title: 'Kinase motifs of the phosphosites: ' + ctxLabel(ctx).toLowerCase(), rows, extraCol: true, heads: ['Kinase (motif class)', 'match / sites', '', 'log2 OR'],
      subHtml: 'Observed phospho-S/T in ' + esc(ctxLabel(ctx).toLowerCase()) + ' against the other phospho-S/T of the same proteins: the share each kinase’s motif matches (Johnson et al. 2023 matrices, percentile ≥ ' + OV.kinase_match +
        '). Within protein (Mantel–Haenszel), BH across the kinases. A prediction from sequence: the residues around a site set its motif, so this says which kinases <i>can</i> write a context, not which do.' + docLink('kinases'),
      subText: 'Phosphosites in ' + ctxLabel(ctx) + ' against the other phosphosites of the same proteins: share matching each kinase motif (percentile >= ' + OV.kinase_match + '), within protein, BH.',
      axis: 'log2 odds ratio, within protein (phosphosites in the context against elsewhere)', valueLabel: 'log2_mh_odds_ratio',
      tsvCols: [['mh_or', (r) => (r.s || r.w).mh_or], ['mh_p', (r) => (r.s || r.w).mh_p], ['motif_class', (r) => r.note]] };
    if (R.sel) spec.selTest = selTestOf(R.sel, 'phosphosites', 'the phosphosites of your selection', 'Mantel–Haenszel, each protein a stratum, BH across the kinases');
    return spec;
  }
  function classSpec() {
    const label = (k) => (OV.classes.find((c) => c[0] === k) || [k, k])[1];
    const what = (k) => (OV.classes.find((c) => c[0] === k) || [k, k, ''])[2];
    const W = tableRows(OV.protein_sets);
    const mk = (pole) => W.filter((r) => r.pole === pole).map((r) => ({ key: pole + ':' + r.class, label: label(r.class), sub: pole, value: l2(r.fold), q: r.q, sig: SIG(r.q),
      extra: fmt(r.k) + '/' + fmt(r.M), w: r, title: '<b>' + esc(label(r.class)) + '</b> (' + esc(what(r.class)) + ')<br>' + fmt(r.k) + ' of the ' + fmt(r.n) + ' proteins with ' + pole + '-pole windows, against ' + fmt(r.M) + ' of ' + fmt(r.N) +
        ' reviewed proteins<br>fold ' + fx(r.fold) + ', one-sided hypergeometric q ' + fmtP(r.q) }));
    let rows;
    const protSel = R.sel && R.sel.prots;   // an evidence-only selection does not change which proteins are counted
    if (protSel) {
      const res = ST.proteinSets(R.prots);
      rows = res.map((x) => ({ key: x.class, label: label(x.class), value: x.k ? l2(x.fold) : null, q: x.q, sig: SIG(x.q), extra: fmt(x.k) + '/' + fmt(x.M), s: x,
        title: '<b>' + esc(label(x.class)) + '</b><br>' + fmt(x.k) + ' of your ' + fmt(x.n) + ' proteins, against ' + fmt(x.M) + ' of ' + fmt(x.N) + ' reviewed proteins<br>fold ' + (x.fold ? fx(x.fold) : '—') + ', q ' + fmtP(x.q) }));
    } else rows = [{ ph: 'Proteins with an acidic-pole window' }].concat(mk('acidic'), [{ ph: 'Proteins with a basic-pole window' }], mk('basic'));
    // section headings ride on the first row of each section, as hbars draws them
    const out = [];
    for (let i = 0; i < rows.length; i++) { if (rows[i].ph) { if (rows[i + 1]) rows[i + 1].head = rows[i].ph; continue; } out.push(rows[i]); }
    for (const r of out) r.main = { or: r.value == null ? null : 2 ** r.value, q: r.q, sig: r.sig, fold: true };
    return { id: 'ch-classes', title: 'Kinds of protein', rows: out, fold: true, extraCol: true, heads: ['Class', 'in set / all', '', 'log2 fold'],
      subHtml: (protSel ? 'Your proteins' : 'Proteins with a window in the top tenth of each pole') + ' against all ' + fmt(B.counts.proteins) + ' reviewed human proteins: one-sided hypergeometric test, fold enrichment (k/n ÷ M/N), BH across the classes.' + docLink('tests'),
      subText: (protSel ? 'Your proteins' : 'Proteins with a pole window') + ' against all reviewed human proteins: one-sided hypergeometric, fold enrichment, BH.',
      axis: 'log2 fold enrichment vs all reviewed human proteins', valueLabel: 'log2_fold_enrichment' };
  }
  const CHARGE_ROWS = ['ph', 'py', 'kac', 'kacyl', 'cit', 'adpr'];
  function chargeSpec() {
    const W = new Map(tableRows(OV.charge_shift).map((r) => [r.kind, r]));
    // the charge changes are computed with every observed site, so only a selection of proteins retests them
    const sel = !!(R.sel && R.sel.prots), prots = sel ? R.sel.prots : null;
    const rows = CHARGE_ROWS.filter((k) => W.has(k)).map((k) => {
      const w = W.get(k), x = sel ? ST.chargeShift(k, prots) : w;
      return { key: k, label: kindLabel(k), value: x.above && x.below ? +Math.log2(x.above / x.below).toFixed(3) : null, p: x.sign_p, w, s: sel ? x : null,
        ref: sel && w.above && w.below ? +Math.log2(w.above / w.below).toFixed(3) : undefined, extra: fmt(x.above) + '/' + fmt(x.above + x.below) };
    });
    C.bhAdjust(rows, (r) => r.p);
    if (sel) { const wq = rows.map((r) => ({ p: r.w.sign_p })); C.bhAdjust(wq); rows.forEach((r, i) => { r.refQ = wq[i].q; }); }
    for (const r of rows) {
      const x = r.s || r.w;
      r.sig = SIG(r.q);
      r.title = '<b>' + esc(r.label) + '</b><br>In ' + fmt(x.above) + ' of ' + fmt(x.above + x.below) + ' IDRs the observed sites raise δ (charge blocks) more than as many ' + esc(r.label.split('-')[0]) +
        ' residues placed at random in the same IDR; in ' + fmt(x.below) + ' less<br>sign test q ' + fmtP(r.q) + ' · δ rises at all in ' + Math.round(100 * (x.share_delta_up || 0)) + '% of these IDRs';
      r.main = { or: r.value == null ? null : 2 ** r.value, q: r.q, sig: r.sig };
    }
    const spec = { id: 'ch-charge', title: 'Do modifications sharpen or soften charge blocks?', rows, extraCol: true, heads: ['Modification', 'sharper / IDRs', '', 'log2 odds'],
      subHtml: 'Each IDR with its observed sites modified (phosphate −2, acetyl-K 0, acidic acyl-K −1, citrulline 0, ADP-ribosyl-S −2; idrkit.ptm), against as many residues of the same type in the same IDR modified at random (20 draws): the odds that the observed sites make the charge blocks sharper (δ, Das & Pappu) than random placement does. Sign test, BH across the modifications.' + docLink('charge'),
      subText: 'Odds that the observed sites raise delta more than random placement of as many residues of the same type in the same IDR; sign test, BH.',
      axis: 'log2 odds: observed sites sharpen charge blocks more than random placement', valueLabel: 'log2_odds_sharper_than_random',
      tsvCols: [['idrs_sharper', (r) => (r.s || r.w).above], ['idrs_softer', (r) => (r.s || r.w).below], ['sign_p', (r) => r.p]] };
    if (sel) spec.selTest = { n: rows.reduce((s, r) => s + (r.s ? r.s.idrs : 0), 0), N: 0, unit: 'IDRs', what: selectionText(), text: 'the modified IDRs of your proteins', method: 'sign test, BH across the modifications' };
    return spec;
  }
  function crosstalkSpec() {
    const rows = tableRows(OV.crosstalk).map((r) => ({ key: r.offset, label: (r.offset > 0 ? '+' : '') + r.offset + ' from a phospho-S/T', value: l2(r.mh_or), q: r.q, sig: SIG(r.q),
      extra: fmt(r.methyl_r_at_offset), title: '<b>Arginines ' + Math.abs(r.offset) + ' residues ' + (r.offset < 0 ? 'before' : 'after') + ' an observed phospho-S/T</b><br>' + fmt(r.methyl_r_at_offset) + ' methylated, ' + fmt(r.other_r_at_offset) +
        ' not<br>within protein and disorder: OR ' + C.fmtOr(r.mh_or, r.mh_lo, r.mh_hi) + ', q ' + fmtP(r.q), main: { or: r.mh_or, lo: r.mh_lo, hi: r.mh_hi, q: r.q, sig: SIG(r.q) } }));
    return { id: 'ch-crosstalk', title: 'Methylarginines around phosphosites', rows, keepOrder: true, extraCol: true, heads: ['Arginine at', 'methyl-R', '', 'log2 OR'],
      subHtml: 'Arginines at each offset from an observed phospho-S/T: methylated against not, within protein and disorder (each protein’s arginines inside IDRs and outside them are strata, as both marks favour IDRs). Mantel–Haenszel, BH across the offsets. Whole dataset only.' + docLink('tests'),
      subText: 'Arginines at each offset from an observed phosphosite, methylated against not, within protein and disorder (Mantel-Haenszel), BH.',
      axis: 'log2 odds ratio, within protein and disorder', valueLabel: 'log2_mh_odds_ratio' };
  }
  function partnersSpec() {
    const rows = [];
    for (const [net, name] of [['string700', 'STRING physical ≥ 700'], ['complexes', 'Complex Portal']]) {
      const sub = tableRows(OV.partners).filter((r) => r.network === net && !r.with_histones);
      sub.forEach((r, i) => rows.push({ key: net + r.pair, head: i ? null : name, sub: net === 'complexes' ? 'Complex Portal' : 'STRING', label: r.pair.replace('contrast: acidic-basic vs same-charge', 'opposite vs same charge'), value: l2(r.obs_over_exp), q: r.p_two_sided, sig: SIG(r.p_two_sided),
        extra: fmt(r.observed), title: '<b>' + esc(name) + ': ' + esc(r.pair) + '</b><br>observed ' + fmt(r.observed) + ', expected ' + (r.expected == null ? '—' : fmt(Math.round(r.expected))) + ' (shuffled within degree and nuclear status)<br>O/E ' + (r.obs_over_exp ?? 0).toFixed(2) + ', p ' + fmtP(r.p_two_sided),
        main: { or: r.obs_over_exp, q: r.p_two_sided, sig: SIG(r.p_two_sided) } }));
    }
    return { id: 'ch-partners', title: 'Do the poles pair across binding partners?', rows, keepOrder: true, extraCol: true, heads: ['Partner pair', 'edges', '', 'log2 O/E'],
      subHtml: 'Partner pairs by the pole class of their proteins (acidic, basic, bipolar, neither), observed against 2,000 shuffles of the classes within degree decile and nuclear status; histone-like proteins left out. The contrast compares opposite-charge pairs with same-charge pairs. Two-sided permutation p. Whole dataset only.' + docLink('tests'),
      subText: 'Partner pairs by pole class, observed over expected under class shuffles within degree and nuclear status; histones left out.',
      axis: 'log2 observed / expected partner pairs', valueLabel: 'log2_obs_over_exp' };
  }
  function stoichSpec(which) {
    const label = which === 'occ' ? 'Phosphosite occupancy (Koenig et al. 2026)' : 'Acetyl-K stoichiometry (Hansen et al. 2019)';
    const measure = which === 'occ' ? 'phospho-S/T occupancy (Koenig 2026)' : 'acetyl-K stoichiometry (Hansen 2019)';
    const ctxs = ['acidic_pole', 'basic_pole', 'window', 'idr', 'blocky', 'basic_patch', 'domain'];
    const W = tableRows(OV.stoichiometry).filter((r) => r.measure === measure);
    const rows = ctxs.map((c) => { const w = W.find((r) => r.context === c); const x = R.sel ? ST.stoich(which, c, R.sel) : w; return { c, w, x }; }).filter((r) => r.w);
    const f = which === 'occ' ? (v) => (100 * v).toFixed(v < 0.1 ? 1 : 0) + '%' : (v) => (100 * v).toFixed(v < 0.001 ? 3 : 2) + '%';
    const box = [];
    for (const r of rows) {
      box.push([ctxLabel(r.c) + ': in', r.x.q_in, 'in', 'Mann–Whitney p ' + fmtP(r.x.mwu_p) + ' against outside']);
      box.push(['outside', r.x.q_out, '', '']);
    }
    return { id: 'ch-' + which, title: label, rows: [], figSvg: () => boxSvg(box, which === 'occ' ? 'occupancy, % of the protein (median of the cell lines, unstimulated)' : 'acetylation stoichiometry, % of the protein (HeLa)', f), figCss: BOX_CSS,
      subHtml: (which === 'occ' ? 'Phospho-S/T occupancy measured in 33 cell lines (median per site)' : 'Acetyl-K stoichiometry measured in HeLa') + ': sites inside each context against sites outside it' + (R.sel ? ' in your proteins' : '') + ', Mann–Whitney test. Pink: inside.' + docLink('occupancy'),
      subText: label + ': inside each context against outside, Mann-Whitney.',
      figTsv: () => '# ' + label + (R.sel ? ' · selection: ' + selectionText() : '') + '\ncontext\tside\tp05\tp25\tmedian\tp75\tp95\tn\tmwu_p\n' + rows.flatMap((r) => [['in', r.x.q_in], ['out', r.x.q_out]].filter(([, q]) => q).map(([s, q]) => [r.c, s].concat(q, [r.x.mwu_p]).join('\t'))).join('\n') + '\n', meta: rows };
  }

  // ---------------------------------------------------------------- is a domain flank a unit? (dtidr_agnostic_profile.py)
  // Every IDR read from each end; each depth against the same IDR's residues shuffled. Series: the ends that touch a
  // DNA-binding domain, an RNA-binding domain, any other domain, and the protein's own terminus.
  const AG_SERIES = [['domain:dna', 'DNA-binding domain (TF)', 's1'], ['domain:rna', 'RNA-binding domain', 's2'], ['domain:other', 'any other domain', 's3'], ['terminus', 'protein terminus (no domain)', 's4']];
  const AG_MEASURES = [['basic', 'K/R fraction', 100, 'percentage points of K or R'], ['acidic', 'D/E fraction', 100, 'percentage points of D or E'],
    ['net_charge', 'net charge', 1, 'net charge per residue'], ['charge_bias', 'charge concentration', 1, '|net charge| of the 10 residues from that depth'],
    ['mod_ph', 'phosphosites', 100, 'percentage points of S/T phosphorylated'], ['mod_mR', 'methylarginines', 100, 'percentage points of R methylated'], ['mod_ub', 'ubiquitin sites', 100, 'percentage points of K ubiquitinated']];
  function agSvg(measure, compact) {
    const [, label, scale, unit] = AG_MEASURES.find((m) => m[0] === measure);
    const rows = tableRows(OV.agnostic).filter((r) => r.measure === measure);
    const bins = [...new Map(rows.map((r) => [r.bin, r.depth_lo])).entries()].sort((x, y) => x[1] - y[1]).map((x) => x[0]);
    const ser = AG_SERIES.map(([key, name, cls]) => ({ key, name, cls, pts: bins.map((b) => rows.find((r) => r.end === key && r.bin === b) || null) }));
    const W = compact ? 300 : 720, H = compact ? 118 : 340, L = compact ? 8 : 64, Rm = compact ? 8 : 190, T = 12, Bm = compact ? 16 : 46;
    const vals = ser.flatMap((x) => x.pts.filter(Boolean).flatMap((p) => (compact ? [p.excess] : [p.ci_low, p.ci_high]))).map((v) => v * scale).concat([0]);
    let lo = Math.min(...vals), hi = Math.max(...vals);
    const pad = (hi - lo) * 0.08 || 1; lo -= pad; hi += pad;
    const X = (i) => L + (bins.length === 1 ? 0 : (i / (bins.length - 1)) * (W - L - Rm)), Y = (v) => T + (1 - (v - lo) / (hi - lo)) * (H - T - Bm);
    let g = '';
    if (!compact) {
      const step = [0.25, 0.5, 1, 2, 2.5, 5, 10].map((x) => (scale === 1 ? x / 100 : x)).find((x) => (hi - lo) / x <= 6) || 10;
      for (let v = Math.ceil(lo / step) * step; v <= hi + 1e-9; v += step) g += '<line class="ag-grid" x1="' + L + '" x2="' + (W - Rm) + '" y1="' + Y(v).toFixed(1) + '" y2="' + Y(v).toFixed(1) + '"/><text class="ag-tick" x="' + (L - 6) + '" y="' + (Y(v) + 4).toFixed(1) + '" text-anchor="end">' + (scale === 1 ? v.toFixed(2) : (v > 0 ? '+' : '') + v.toFixed(Math.abs(step) < 1 ? 1 : 0)) + '</text>';
      bins.forEach((b, i) => { g += '<text class="ag-tick" x="' + X(i).toFixed(1) + '" y="' + (H - Bm + 16) + '" text-anchor="middle">' + b + '</text>'; });
      g += '<text class="ag-axt" x="' + ((L + W - Rm) / 2).toFixed(1) + '" y="' + (H - 6) + '" text-anchor="middle">residues from the IDR’s end (1 touches the domain or terminus)</text>';
      g += '<text class="ag-axt" transform="translate(14 ' + ((T + H - Bm) / 2).toFixed(1) + ') rotate(-90)" text-anchor="middle">excess over the same IDR shuffled</text>';
    }
    g += '<line class="ag-zero" x1="' + L + '" x2="' + (W - Rm) + '" y1="' + Y(0).toFixed(1) + '" y2="' + Y(0).toFixed(1) + '"/>';
    const labs = [];
    for (const sr of ser) {
      const pts = sr.pts.map((p, i) => (p ? [X(i), Y(p.excess * scale), p] : null)).filter(Boolean);
      if (!pts.length) continue;
      g += '<polyline class="ag-line ag-' + sr.cls + '" points="' + pts.map(([x, y]) => x.toFixed(1) + ',' + y.toFixed(1)).join(' ') + '"/>';
      for (const [x, y, p] of pts) {
        const sig = p.ci_low > 0 || p.ci_high < 0;
        if (!compact) g += '<line class="ag-ci ag-' + sr.cls + '" x1="' + x.toFixed(1) + '" x2="' + x.toFixed(1) + '" y1="' + Y(p.ci_low * scale).toFixed(1) + '" y2="' + Y(p.ci_high * scale).toFixed(1) + '"/>';
        g += '<circle class="ag-pt ag-' + sr.cls + (sig ? '' : ' ns') + '" cx="' + x.toFixed(1) + '" cy="' + y.toFixed(1) + '" r="' + (compact ? 2.6 : 4) + '" data-tip="' +
          esc('<b>' + esc(sr.name) + ', residues ' + p.bin + '</b><br>' + esc(label) + ' ' + (p.excess * scale >= 0 ? '+' : '') + (p.excess * scale).toFixed(scale === 1 ? 3 : 2) + ' ' + (scale === 1 ? '' : 'points ') + 'over the same IDR shuffled<br>95% CI ' +
            (p.ci_low * scale).toFixed(scale === 1 ? 3 : 2) + ' to ' + (p.ci_high * scale).toFixed(scale === 1 ? 3 : 2) + ' (bootstrap over proteins) · ' + fmt(p.ends) + ' IDR ends on ' + fmt(p.proteins) + ' proteins') + '"/>';
      }
      if (!compact) labs.push([pts[pts.length - 1][1], sr]);
    }
    // direct labels at the right, spread so they do not collide
    labs.sort((x, y) => x[0] - y[0]);
    let last = -1e9;
    for (const [y0, sr] of labs) { const y = Math.max(y0, last + 13); last = y; g += '<text class="ag-lab ag-' + sr.cls + '-t" x="' + (W - Rm + 8) + '" y="' + (y + 4).toFixed(1) + '">' + esc(sr.name) + '</text>'; }
    return '<svg class="ag-plot' + (compact ? ' compact' : '') + '" viewBox="0 0 ' + W + ' ' + H + '" role="img" aria-label="' + esc(label + ' by distance from the IDR end, against the same IDR shuffled') + '">' + g + '</svg>';
  }
  const AG_CSS = (pal) => '.ag-grid{stroke:' + pal.grid + '}.ag-zero{stroke:' + pal.muted + ';stroke-width:1}.ag-tick,.ag-axt{font-size:11px;fill:' + pal.muted + '}.ag-lab{font-size:11px}.ag-line{fill:none;stroke-width:2}.ag-ci{stroke-width:1.2;opacity:.6}' +
    '.ag-pt{stroke:#fff;stroke-width:1.5}.ag-pt.ns{fill:#fff!important}' + [['s1', '#2a78d6'], ['s2', '#eb6834'], ['s3', '#1baf7a'], ['s4', '#c98500']].map(([k, c]) => '.ag-' + k + '{stroke:' + c + ';fill:' + c + '}.ag-line.ag-' + k + '{fill:none}.ag-' + k + '-t{fill:' + c + '}').join('');
  function agSpec(measure) {
    const [, label, scale, unit] = AG_MEASURES.find((m) => m[0] === measure);
    const rows = tableRows(OV.agnostic).filter((r) => r.measure === measure && AG_SERIES.some((x) => x[0] === r.end));
    return { id: 'ch-agnostic-' + measure, title: 'Is the sequence beside a domain different? ' + label, rows: [], figSvg: () => agSvg(measure, false), figCss: AG_CSS,
      subHtml: 'Every IDR read from each end, each depth against the same IDR’s residues in random order (20 shuffles), so composition and length are kept and only position is removed. Lines: the IDR ends that touch each kind of domain, and the ends at the protein’s own terminus. Filled points: 95% CI (bootstrap over proteins) excludes 0. Whole dataset only.' + docLink('agnostic'),
      subText: label + ' (' + unit + ') by distance from the IDR end, against the same IDR shuffled; 95% CI by bootstrap over proteins.',
      figTsv: () => '# ' + label + ', ' + unit + '; excess over the same IDR shuffled\nend\tbin\texcess\tci_low\tci_high\tends\tproteins\n' + rows.map((r) => [r.end, r.bin, r.excess * scale, r.ci_low * scale, r.ci_high * scale, r.ends, r.proteins].join('\t')).join('\n') + '\n' };
  }

  // ---------------------------------------------------------------- Overview
  // the tiles: key, group, context or builder, question
  const TILES = [
    ['agnostic', 'The unit', 'Is the sequence beside a folded domain different from the rest of its IDR?'],
    ['acidic_pole', 'The two poles', 'What does the acidic pole carry?'], ['basic_pole', 'The two poles', 'What does the basic pole carry?'], ['blocky', 'The two poles', 'What sits in IDRs with charge blocks?'],
    ['acidic_stretch', 'Sequence elements', 'What sits in acidic stretches?'], ['basic_patch', 'Sequence elements', 'What sits in basic patches?'], ['rg_cluster', 'Sequence elements', 'What sits in RG/RGG clusters?'],
    ['basophilic', 'Sequence elements', 'Are the arginines and lysines just before a phosphosite modified?'],
    ['idr', 'Order and disorder', 'What sits in IDRs?'], ['domain', 'Order and disorder', 'What sits in folded domains?'], ['window', 'Order and disorder', 'What sits beside a domain?'],
    ['kin', 'Writers', 'Which kinases’ motifs do each pole’s phosphosites match?'], ['charge', 'Charge', 'Do modifications sharpen or soften charge blocks?'],
    ['occupancy', 'Occupancy', 'Are sites in the poles more fully modified?'], ['classes', 'Proteins', 'What kinds of protein carry the poles?'],
    ['crosstalk', 'Crosstalk and partners', 'Are arginines next to phosphosites methylated?'], ['partners', 'Crosstalk and partners', 'Do the poles pair across binding partners?'],
  ];
  const KIN_CTX = ['acidic_pole', 'basic_pole', 'blocky', 'acidic_stretch', 'basic_patch', 'idr', 'domain'];
  const tileRows = (rows, n = 6) => {
    const ok = rows.filter((r) => r.value != null && r.main && r.main.or);
    const up = ok.filter((r) => r.sig && r.value > 0).sort((a, b) => b.value - a.value), dn = ok.filter((r) => r.sig && r.value < 0).sort((a, b) => a.value - b.value);
    const nUp = Math.min(up.length, dn.length ? n - 2 : n), nDn = Math.min(dn.length, n - nUp);
    return up.slice(0, nUp).concat(dn.slice(0, nDn)).map((r) => ({ label: r.label, sub: r.sub, main: r.main }));
  };
  const answer = (rows, fold) => {
    const ok = rows.filter((r) => r.value != null && r.sig && Math.abs(r.value) >= Math.log2(C.minEffect()));
    const up = ok.filter((r) => r.value > 0).sort((a, b) => b.value - a.value), dn = ok.filter((r) => r.value < 0).sort((a, b) => a.value - b.value);
    const eff = (r) => (fold ? fx(2 ** r.value) + '-fold' : fx(2 ** r.value) + '×');
    const one = (list) => list.slice(0, 2).map((r) => esc(r.label) + ' <span class="an">' + eff(r) + '</span>').join(', ');
    if (up.length) return C.am(1) + one(up) + (dn.length ? '; ' + C.am(-1) + one(dn) : '');
    if (dn.length) return C.am(-1) + one(dn);
    return C.am(0) + 'none differs at q ≤ ' + PADJ;
  };
  const wholeOf = (spec) => spec.rows.map((r) => ({ ...r, value: spec.selTest ? r.ref : r.value, sig: spec.selTest ? SIG(r.refQ) : r.sig }));
  const SECTIONS = {};
  function buildSection(key) {
    if (CTX_ORDER().includes(key)) {
      const spec = ctxSpec(key);
      return { spec, rows: spec.rows };
    }
    if (key === 'kin') {
      const a = kinaseSpec('acidic_pole'), b = kinaseSpec('basic_pole');
      return { spec: kinaseSpec(UI.kinCtx || 'acidic_pole'), groups: [['Acidic pole', a.rows], ['Basic pole', b.rows]], rows: a.rows.concat(b.rows), specs: [a, b] };
    }
    if (key === 'charge') { const spec = chargeSpec(); return { spec, rows: spec.rows }; }
    if (key === 'classes') { const spec = classSpec(); return { spec, rows: spec.rows, fold: true }; }
    if (key === 'crosstalk') { const spec = crosstalkSpec(); return { spec, rows: spec.rows, whole: true }; }
    if (key === 'partners') { const spec = partnersSpec(); return { spec, rows: spec.rows, whole: true }; }
    if (key === 'occupancy') { const occ = stoichSpec('occ'), ac = stoichSpec('ac'); return { spec: occ, spec2: ac, occupancy: true }; }
    if (key === 'agnostic') return { spec: agSpec(UI.agMeasure || 'basic'), agnostic: true, whole: true };
    return null;
  }
  function tileOf(key, sec) {
    if (sec.agnostic) {
      const r = tableRows(OV.agnostic), get = (end, bin) => r.find((x) => x.measure === 'basic' && x.end === end && x.bin === bin) || {};
      const dna = get('domain:dna', '1-5'), dna2 = get('domain:dna', '16-20'), oth = get('domain:other', '1-5');
      return { vis: agSvg('basic', true) + '<div class="ag-key">' + AG_SERIES.map(([, n, c]) => '<span><i class="ag-sw ag-' + c + '"></i>' + esc(n) + '</span>').join('') + '</div>',
        head: 'Only beside DNA-binding domains: <b>+' + (100 * dna.excess).toFixed(1) + '</b> points of K/R at residues 1–5, gone by 16–20 (' + (100 * dna2.excess).toFixed(1) + '); other domains +' + (100 * oth.excess).toFixed(1) + ' at 1–5 only' };
    }
    if (sec.occupancy) {
      const rows = sec.spec.meta, mx = Math.max(...rows.map((r) => (r.x.median_in || 0)), ...rows.map((r) => r.x.median_out || 0));
      const meters = rows.filter((r) => ['acidic_pole', 'basic_pole', 'idr', 'domain'].includes(r.c)).map((r) => ({ label: ctxLabel(r.c), v: (r.x.median_in || 0) / mx, ref: (r.x.median_out || 0) / mx,
        text: r.x.median_in == null ? '—' : (100 * r.x.median_in).toFixed(1) + '%', tip: '<b>' + esc(ctxLabel(r.c)) + '</b><br>median occupancy ' + (100 * (r.x.median_in || 0)).toFixed(1) + '% inside (' + fmt(r.x.n_in) + ' sites), ' +
          (100 * (r.x.median_out || 0)).toFixed(1) + '% outside (tick); Mann–Whitney p ' + fmtP(r.x.mwu_p) }));
      const idr = rows.find((r) => r.c === 'idr');
      return { vis: C.aMeters(meters, { tight: true }), head: 'Median phospho-S/T occupancy inside, tick outside' + (idr ? '; IDRs ' + (100 * idr.x.median_in).toFixed(1) + '% vs ' + (100 * idr.x.median_out).toFixed(1) + '%' : '') };
    }
    const vis = sec.groups ? C.aForest(sec.groups.map(([g, rows]) => [g, tileRows(rows, 3)]), { compact: true }) : C.aForest([[null, tileRows(sec.rows)]], { compact: true });
    return { vis, head: answer(sec.rows, sec.fold), whole: R.scoped && !sec.whole ? answer(wholeOf(sec.spec), sec.fold) : '' };
  }

  function selPanel() {
    const checks = (name, items, on) => items.map(([k, l]) => '<label class="check"><input type="checkbox" name="' + name + '" value="' + k + '"' + (on.includes(k) ? ' checked' : '') + '>' + esc(l) + '</label>').join('');
    const src = F.src == null ? allSrc() : F.src;
    return '<details class="panel sel-panel"' + (R.scoped ? ' open' : '') + '><summary><b>Your selection</b> <span class="muted">' + esc(selectionText()) + '</span>' +
      (R.scoped ? ' · <b>' + fmt(R.sites.length) + '</b> sites on <b>' + fmt(R.prots.length) + '</b> proteins' : ' · every test below runs again on a selection') + '</summary>' +
      '<form class="filters sel-form" id="selform"><label class="wide">Genes or accessions<input type="text" name="genes" value="' + esc(F.genes) + '" placeholder="e.g. NPM1 NCL HNRNPU, or one prefix"></label>' +
      '<fieldset><legend>Proteins with (all of)</legend>' + checks('has', OV.has.map((h) => [h, HAS_LABEL[h]]), F.has) + '</fieldset>' +
      '<fieldset><legend>Protein class (any of)</legend>' + checks('cls', OV.classes.map((c) => [c[0], c[1]]), F.cls) + '</fieldset>' +
      '<fieldset><legend>Sites from (any of)</legend>' + OV.sources.map(([k, b, l]) => '<label class="check"><input type="checkbox" name="src" value="' + b + '"' + (src & b ? ' checked' : '') + '>' + esc(l) + '</label>').join('') +
      '<label>PTMeXchange-only sites<select name="conf">' + ['any grade', 'Bronze or better', 'Silver or better', 'Gold only'].map((l, k) => '<option value="' + k + '"' + (F.conf === k ? ' selected' : '') + '>' + l + '</option>').join('') + '</select></label></fieldset>' +
      '<div class="sel-actions"><button class="btn" type="submit">Apply</button>' + (R.scoped ? '<button class="btn" type="button" id="selclear">Clear the selection</button>' : '') + '</div></form></details>';
  }
  function bindSelPanel(path) {
    const form = $('#selform');
    if (!form) return;
    const apply = () => {
      const fd = new FormData(form);
      F.genes = (fd.get('genes') || '').trim();
      F.has = fd.getAll('has'); F.cls = fd.getAll('cls');
      const src = fd.getAll('src').reduce((b, v) => b | +v, 0);
      F.src = src === allSrc() ? null : src;
      F.conf = +(fd.get('conf') || 0);
      setParams(path, selParams());
      route();
    };
    form.addEventListener('submit', (e) => { e.preventDefault(); apply(); });
    form.addEventListener('change', (e) => { if (e.target.type === 'checkbox' || e.target.tagName === 'SELECT') apply(); });
    const clear = $('#selclear');
    if (clear) clear.addEventListener('click', () => { Object.assign(F, { genes: '', has: [], cls: [], src: null, conf: 0 }); setParams(path, selParams()); route(); });
  }

  async function overviewView(params, token, sub) {
    if (constraints().length) view().innerHTML = '<p class="loading">Testing your selection…</p>';
    await computeSelection();
    if (token !== routeToken) return;
    UI.charts = {};
    if (sub) { sectionView(sub); return; }
    C.A_ROWS.length = 0;
    const ov = {};
    for (const [key, group, q] of TILES) {
      const sec = buildSection(key);
      if (!sec) continue;
      SECTIONS[key] = sec;
      const title = CTX_ORDER().includes(key) ? ctxLabel(key) : { agnostic: 'Domain flanks, read without a window', kin: 'Kinase motifs', charge: 'Charge blocks', occupancy: 'Occupancy', classes: 'Kinds of protein', crosstalk: 'Crosstalk', partners: 'Partners' }[key];
      const how = sec.whole ? 'whole' : key === 'occupancy' ? 'retest' : 'retest';
      C.foldHtml(key, title, '', false, '', null, { q, how });
      C.A_ROWS[C.A_ROWS.length - 1].group = group;
      ov[key] = tileOf(key, sec);
      if (sec.whole && R.scoped) ov[key].head = 'Whole dataset only · ' + ov[key].head;
    }
    const c = B.counts, rt = tableRows(OV.residue_tests);
    const get = (k, ctx) => rt.find((r) => r.kind === k && r.context === ctx) || {};
    const cs = tableRows(OV.charge_shift);
    const hi = (k, ctx, go, glyph, label) => { const r = get(k, ctx); return { go, glyph, big: fx(r.mh_or) + '×', label, kind: kindLabel(k) + ' · within protein', dir: r.mh_or >= 1 ? 'up' : 'down',
      tip: '<b>' + esc(kindLabel(k)) + ' in ' + esc(ctxLabel(ctx).toLowerCase()) + '</b><br>' + fmt(r.sites_in) + ' of ' + fmt(r.sites) + ' sites; OR ' + C.fmtOr(r.mh_or, r.mh_lo, r.mh_hi) + ' within protein, q ' + fmtP(r.q_in_context) }; };
    const ph = cs.find((r) => r.kind === 'ph') || {}, ad = cs.find((r) => r.kind === 'adpr') || {};
    const lead = R.scoped ? { q: 'Your selection', pic: '<p class="a-sel-what">' + esc(selectionText()) + ' · <b>' + fmt(R.sites.length) + '</b> sites on <b>' + fmt(R.prots.length) + '</b> proteins</p>',
      verdict: R.prots.length < 10 ? '<span class="a-note">Too few proteins for the tests to mean much (fewer than 10); the numbers are shown as they come.</span>' : '' }
      : { q: 'Where do the modifications that change charge sit?', hi: C.aHighlights([
        hi('ck2', 'acidic_stretch', 'acidic_stretch', C.aGlyph.motif('SEEDEE'), 'CK2 sites in acidic stretches'),
        hi('mR', 'rg_cluster', 'rg_cluster', C.aGlyph.motif('RGGRGG'), 'Methylarginines in RG/RGG clusters'),
        hi('cit', 'rg_cluster', 'rg_cluster', C.aGlyph.motif('RGGRGG'), 'Citrullines in RG/RGG clusters'),
        hi('adpr', 'basic_patch', 'basic_patch', C.aGlyph.motif('KSKKSK'), 'ADP-ribosylserines in basic patches'),
        { go: 'charge', glyph: C.aGlyph.pins(), big: Math.round(100 * ph.above / (ph.above + ph.below)) + '%', label: 'IDRs where phosphorylation sharpens charge blocks more than random sites would', kind: 'sign test, p ' + fmtP(ph.sign_p), dir: 'up',
          tip: 'Phosphorylating the observed sites raises δ more than phosphorylating as many random S/T in the same IDR in ' + fmt(ph.above) + ' of ' + fmt(ph.above + ph.below) + ' IDRs; ADP-ribosylation the reverse (' + fmt(ad.below) + ' of ' + fmt(ad.above + ad.below) + ' softer)' },
      ]) + '<p class="a-key a-key-line"><span><i class="aq-pt up"></i>more often</span><span><i class="aq-pt down"></i>less often</span><span><i class="aq-pt same"></i>about the same (within ' + C.minEffect() + '×)</span>' +
        '<span class="muted">In each tile: modified residues against unmodified residues of the same type in the same proteins, within protein; each opens its chart, numbers and method.</span></p>' };
    view().innerHTML = selPanel() + '<div class="stats">' + [
      [c.windows, 'dtIDR windows (' + fmt(c.acidic_pole_windows) + ' acidic, ' + fmt(c.basic_pole_windows) + ' basic top tenth)'], [c.idrs, 'IDRs (≥ 20 aa)'], [c.sites, 'modification sites'],
      [c.by_kind.ph, 'phospho-S/T'], [c.by_kind.mR + c.by_kind.cit, 'methylated or citrullinated R'], [c.by_kind.kac + c.by_kind.ub + c.by_kind.sumo, 'acetyl, ubiquitin, SUMO-K'], [c.by_kind.adpr + c.by_kind.glu, 'ADP-ribosyl and glutamyl sites'],
    ].map(([n, l]) => '<div class="stat"><b>' + fmt(n) + '</b><span>' + esc(l) + '</span></div>').join('') + '</div>' + C.aHome(lead, ov) + ideasHtml();
    $$('.a-tile[data-go], .a-hi-c[data-go]').forEach((b) => b.addEventListener('click', () => { const s = selParams().toString(); location.hash = 'overview/' + b.dataset.go + (s ? '?' + s : ''); }));
    bindSelPanel('overview');
  }

  function sectionView(key) {
    const sec = buildSection(key);
    if (!sec) { view().innerHTML = '<div class="panel">No section called ' + esc(key) + '. <a href="#overview">Back to the Overview</a></div>'; return; }
    const idx = TILES.findIndex((t) => t[0] === key), t = TILES[idx];
    const s = selParams().toString(), q = s ? '?' + s : '';
    const nav = '<nav class="sec-nav"><a href="#overview' + q + '">‹ All tiles</a>' + (idx > 0 ? '<a href="#overview/' + TILES[idx - 1][0] + q + '">‹ ' + esc(TILES[idx - 1][2]) + '</a>' : '') +
      (idx < TILES.length - 1 ? '<a href="#overview/' + TILES[idx + 1][0] + q + '">' + esc(TILES[idx + 1][2]) + ' ›</a>' : '') + '</nav>';
    let body = '';
    if (key === 'kin') {
      const ctx = UI.kinCtx || 'acidic_pole';
      sec.spec.tools = '<label>Context <select id="kinctx">' + KIN_CTX.map((c) => '<option value="' + c + '"' + (c === ctx ? ' selected' : '') + '>' + esc(ctxLabel(c)) + '</option>').join('') + '</select></label>';
      body = card(sec.spec);
    } else if (sec.occupancy) body = card(sec.spec) + card(sec.spec2);
    else if (sec.agnostic) {
      const m = UI.agMeasure || 'basic';
      sec.spec.tools = '<label>Measure <select id="agm">' + AG_MEASURES.map(([k, l]) => '<option value="' + k + '"' + (k === m ? ' selected' : '') + '>' + esc(l) + '</option>').join('') + '</select></label>';
      body = card(sec.spec, '<div class="ag-key">' + AG_SERIES.map(([, n, c]) => '<span><i class="ag-sw ag-' + c + '"></i>' + esc(n) + '</span>').join('') + '</div>') + agSummaryHtml();
    }
    else body = card(sec.spec, CTX_ORDER().includes(key) ? '<p class="note">Each modification’s residues: ' + OV.kinds.map((k) => esc(k[1]) + ' ' + k[2].split('').join('/')).join(', ') + '. ADP-ribosylation is tested on serines; ubiquitin and SUMO on lysines.</p>' : '');
    view().innerHTML = selPanel() + nav + '<h2 class="sec-q">' + esc(t[2]) + '</h2>' + body + sectionNotes(key) + nav;
    $$('.hbars').forEach((el) => { const id = el.closest('.chart-card').id, spec = UI.charts[id]; if (spec) C.bindBars(el, spec.rows, () => {}); });
    bindExports(view());
    const am = $('#agm');
    if (am) am.addEventListener('change', () => { UI.agMeasure = am.value; sectionView('agnostic'); });
    const kc = $('#kinctx');
    if (kc) kc.addEventListener('change', () => { UI.kinCtx = kc.value; sectionView('kin'); });
    bindSelPanel('overview/' + key);
  }
  function agSummaryHtml() {
    const S = tableRows(OV.agnostic_summary), cp = S.find((r) => r.end === 'domain change points');
    const ends = (k) => (tableRows(OV.agnostic).find((r) => r.end === k && r.measure === 'basic' && r.bin === '1-5') || {}).ends;
    return '<div class="panel note-panel"><p><b>What this shows.</b> Read without drawing a window, the residues beside most folded domains are not a distinct segment: beyond the first five residues (a small basic cap) they look like the rest of their IDR, and protein termini show the same small effects. ' +
      'The exception is the DNA-binding domain of a transcription factor: its ends (' + fmt(ends('domain:dna')) + ') carry a basic extension that fades over about 15 residues, the basic tails known to tune DNA binding. Histone-chaperone folds are too few to read here (' + fmt(ends('domain:chaperone')) + ' ends), with acidic sequence 30–75 residues out. ' +
      (cp ? 'Splitting each domain-touching IDR where its charge changes most finds a split stronger than in shuffles in ' + Math.round(100 * cp.excess_1_5) + '% of ' + fmt(cp.ends) + ' ends (5% expected by chance), at a median of ' + cp.half_depth + ' residues. ' : '') +
      'So a dtIDR window is a coordinate frame (an address next to a fold) whose content varies by domain class, not a universal element; the claims on this site are tests within that frame.</p></div>';
  }
  function sectionNotes(key) {
    const N = {
      acidic_pole: 'The acidic pole is defined from sequence alone (disorder, acidic charge, glutamate tracts); no modification enters the score. CK2 and glutamylation are the expected marks; read methyl-R and acetyl-K as the basic pole’s marks seen from the other side.',
      basic_pole: 'The basic pole is the MicroDomains dtIDR score (disorder, net positive charge, local R/K clustering, RG/GR density); no modification enters it.',
      charge: 'A modification that adds negative charge next to an acidic block lengthens it and raises δ; one inside a basic patch interrupts it and lowers δ. Random placement of as many residues of the same type in the same IDR is the null, so a composition change alone scores zero. Glutamylation is left out: it deepens an existing negative charge, which δ (counted by sign) does not see. Ubiquitin and SUMO attach a protein and are left out.',
      kin: 'Motif matches are predictions: a residue whose neighbours are acidic matches acidophilic kinases whether or not they act on it. The test therefore compares observed phosphosites with other observed phosphosites of the same proteins; it says which kinases could write a context’s sites.',
      crosstalk: 'Methylarginines are enriched at −2 and +2 to +5 from phosphosites, not at the −3 position the basophilic R-x-x-S/T motif needs: shared sequence context (S/G/R-rich segments) more than a kinase-motif link.',
      partners: 'Both poles are enriched in nuclear and chromatin proteins, which bind each other, so the null keeps degree and nuclear status. Like binds like: basic–basic pairs are over-represented, acidic–basic pairs are not.',
      occupancy: 'Occupancy is the share of a site’s protein copies that carry the modification. Mann–Whitney compares the distributions; read the medians with the counts, as the poles hold few measured sites.',
      classes: 'Classes come from the windows table of MicroDomains (TFs: Lambert et al. 2018; RNA-binding domains; curated histone chaperones; histone folds; HMG proteins).',
    };
    return N[key] ? '<div class="panel note-panel"><p>' + esc(N[key]) + '</p></div>' : '';
  }

  // The ideas the numbers support, each pointing to its tile; numbers from overview.json.
  function ideasHtml() {
    const rt = tableRows(OV.residue_tests), g = (k, c) => rt.find((r) => r.kind === k && r.context === c) || {};
    const or = (k, c) => fx(g(k, c).mh_or) + '×';
    const go = (k, l) => '<a href="#overview/' + k + '">' + esc(l) + ' →</a>';
    const ideas = [
      ['neg', 'Tested: the unit', 'Beside most domains there is no distinct segment', 'Read without a window, the residues next to a fold resemble the rest of their IDR except for a five-residue basic cap; termini show the same. Only DNA-binding domains carry a basic extension, fading over about 15 residues. The dtIDR window is therefore an address (a coordinate frame beside a fold), not a universal element.', [go('agnostic', 'Domain flanks, read without a window')]],
      ['data', 'The frame', 'Scored from sequence, tested against modifications', 'Every domain edge of the reviewed proteome is scored from sequence alone, once for basic and once for acidic character; no modification enters either score, so each tile is a test of the scores against data they never saw.', [go('acidic_pole', 'Acidic pole'), go('basic_pole', 'Basic pole')]],
      ['new', 'Two chemistries', 'Each pole has its own writers', 'CK2 sites sit at the acidic pole (' + or('ck2', 'acidic_pole') + ') and in acidic stretches (' + or('ck2', 'acidic_stretch') + '); methylarginines in RG/RGG clusters (' + or('mR', 'rg_cluster') + ') and at the basic pole (' + or('mR', 'basic_pole') + '). The phosphosites of each pole match the motifs of different kinases: acidophilic (CK2, GRK, FAM20C) at one, basophilic (PKC, AKT, PKA) at the other.', [go('kin', 'Kinase motifs')]],
      ['new', 'Charge tuning', 'Modifications are placed to change charge blocks', 'Phosphorylation and citrullination at their observed sites sharpen charge blocks more than the same number of random residues would; ADP-ribosylation softens them, as its −2 lands inside basic patches (' + or('adpr', 'basic_patch') + '). Placement, not just composition, is under selection.', [go('charge', 'Charge blocks')]],
      ['new', 'Arginine', 'PADs and PRMTs pick different arginines', 'Citrullines avoid RG/RGG clusters (' + or('cit', 'rg_cluster') + ') where methylarginines crowd (' + or('mR', 'rg_cluster') + '); both favour IDRs. The two arginine marks compete for few residues.', [go('rg_cluster', 'RG/RGG clusters')]],
      ['neg', 'Tested: not a rule', 'The poles do not pair across partners', 'Opposite-charge partner pairs are not over-represented after degree and nuclear status are held; same-charge pairs are. Charge pairing is real in named systems (chaperone and histone), not who binds whom across the proteome.', [go('partners', 'Partners')]],
      ['hyp', 'Hypothesis', 'A dtIDR is an avidity tuner', 'The folded domain chooses the site; the adjacent residues make a second, low-specificity contact with the same ligand surface. A modification placed to sharpen or soften a charge block retunes that contact, so it should change residence time rather than site selection.', [go('charge', 'Charge blocks')]],
      ['open', 'Open', 'What is not settled', 'Whether CK2 writes the acidic IDRs it matches or only marks them; why phosphosites in the acidic pole are not more occupied; a proteome-wide glutamylation site table (none exists); FINCHES predictions of IDR–IDR attraction for the partner test.', []],
    ];
    return '<section class="ideas-sec"><h2>Ideas</h2><div class="ideas">' + ideas.map(([tag, tl, h, p, links]) => '<article class="idea"><span class="tag ' + tag + '">' + esc(tl) + '</span><h3>' + esc(h) + '</h3><p>' + p + '</p>' + (links.length ? '<div class="foot">' + links.join('') + '</div>' : '') + '</article>').join('') + '</div></section>';
  }

  // ---------------------------------------------------------------- Definitions
  async function definitionsView() {
    const ctx = OV.contexts.map(([k, l, d]) => '<dt id="def-ctx-' + k + '">' + esc(l) + '</dt><dd>' + esc(d) + '.</dd>').join('');
    const kinds = OV.kinds.map(([k, l, res, mod]) => '<tr><td>' + esc(l) + '</td><td>' + res.split('').join(', ') + '</td><td>' + esc(mod || '—') + '</td><td class="num">' + fmt(B.counts.by_kind[k] || 0) + '</td></tr>').join('');
    const S = OV.settings;
    view().innerHTML = '<div class="about defs"><h2>Definitions</h2>' +
      '<div class="panel"><h3 id="def-unit">The unit</h3><dl>' +
      '<dt>dtIDR window</dt><dd>A domain-tethered IDR: the residues beside one edge of a folded domain (InterPro domain or repeat), up to 50 on that side or fewer where the next domain begins. ' + fmt(B.counts.windows) + ' windows on ' + fmt(B.counts.proteins_with_windows) + ' proteins.</dd>' +
      '<dt>Basic score</dt><dd>MicroDomains’ dtIDR score (<code>dtidr_score_z</code>): disorder, net positive charge, local R/K clustering and RG/GR density, each as a percentile, averaged and standardised. Shown as a z-score.</dd>' +
      '<dt>Acidic score</dt><dd>The same construction with the charge inverted: disorder, acidic charge fraction and density of glutamates in acidic runs (<code>acidic_dtidr_score</code>).</dd>' +
      '<dt id="def-poles">Pole</dt><dd>A window in the top tenth of all windows (' + (100 * S.pole_top).toFixed(0) + 'th percentile) for a score: ' + fmt(B.counts.acidic_pole_windows) + ' acidic-pole and ' + fmt(B.counts.basic_pole_windows) + ' basic-pole windows. A window can be both (bipolar).</dd>' +
      '<dt>IDR</dt><dd>A MobiDB disordered region (overlapping fragments merged) of ' + S.min_idr + ' or more residues: ' + fmt(B.counts.idrs) + '.</dd></dl></div>' +
      '<div class="panel"><h3 id="def-agnostic">Reading a flank without a window</h3><p>Every IDR is read residue by residue from each of its ends. At each depth, a measure (K/R fraction, D/E fraction, net charge, the |net charge| of the next 10 residues, or the share of residues of a type that carry a modification) is compared with the same IDR’s residues in random order (20 shuffles), so length and composition are kept and only position is removed. Ends are classed by what lies beyond them: a folded domain (InterPro domain or repeat within 10 residues; split by the domain’s class), the protein’s terminus, or ordered sequence that is not a domain. 95% intervals: bootstrap over proteins. A per-end change point (the split with the largest proximal–distal difference in charge concentration, against shuffles) gives the length scale without choosing one.</p></div>' +
      '<div class="panel"><h3 id="def-contexts">Where a residue sits</h3><p>Each test asks whether a modification’s residues sit in a context more often than the unmodified residues of the same type in the same proteins. A residue can be in several contexts; each is tested on its own, inside against outside.</p><dl>' + ctx + '</dl></div>' +
      '<div class="panel"><h3 id="def-kinds">Modifications</h3><p>Sites as the sources report them, placed on the UniProt 2026_03 reviewed sequences by residue and, where the source gives peptides, by peptide. The charge model is the one the charge tests use (idrkit.ptm).</p>' +
      '<div class="table-wrapper" tabindex="0" role="region" aria-label="Modification kinds"><table class="data"><thead><tr><th>Modification</th><th>Residues tested</th><th>Charge model</th><th class="num">Sites</th></tr></thead><tbody>' + kinds + '</tbody></table></div>' +
      '<p class="note">Charge of a modified residue at pH 7.4: phospho −2; acetyl-K and neutral acyls 0; succinyl, malonyl and glutaryl-K −1; citrulline 0; ADP-ribosyl-S −2 (on D or E: −2 in place of −1); one added glutamate −2 in place of −1 (a floor); methylation unchanged. Ubiquitin and SUMO are left out of charge.</p></div>' +
      '<div class="panel"><h3 id="def-charge">Charge measures</h3><dl>' +
      '<dt>NCPR, FCR</dt><dd>Net charge per residue and fraction of charged residues (K, R +1; D, E −1; H 0).</dd>' +
      '<dt>δ (delta)</dt><dd>Das & Pappu’s charge patterning (PNAS 110:13392, 2013): the mean over all blobs of 5 and of 6 residues of (σ − σ<sub>blob</sub>)², σ = NCPR²/FCR, as localCIDER computes it (idrkit.delta). Higher: charges sit in local blocks.</dd>' +
      '<dt>δ z</dt><dd>δ against 40 shuffles of the same residues (length and composition kept). Above 0: blockier than the composition gives. An IDR is <b>blocky</b> at δ z ≥ ' + S.blocky_z + '.</dd>' +
      '<dt>SCD z</dt><dd>Sequence charge decoration (Sawle & Ghosh 2015) against the same shuffles: long-range segregation; below 0, more segregated.</dd>' +
      '<dt>Sharper or softer than random</dt><dd>For each IDR and modification, δ with the observed sites modified against δ with as many residues of the same type in the same IDR modified at random (20 draws). The share of IDRs where the observed change is larger is tested with a sign test.</dd></dl></div>' +
      '<div class="panel"><h3 id="def-kinases">Kinase motifs</h3><p>Every serine and threonine in an IDR or window, and every observed phosphosite, scored with the Ser/Thr kinase matrices of Johnson et al. 2023 (Nature 613:759) through kinase-library 1.8.0: for each kinase, the percentile of the residue’s score among the reference human phosphoproteome. A <b>motif match</b> is a percentile of ' + OV.kinase_match + ' or more. Representatives name each motif class: ' +
      OV.kinase_rows.map((r) => esc(r[1]) + ' (' + esc(r[2]) + ')').join(', ') + '. Strong unobserved CK2 motifs (CK2 percentile ≥ 99) are drawn on the protein maps: ' + fmt(B.counts.dark_ck2_motifs) + '.</p></div>' +
      '<div class="panel"><h3 id="def-occupancy">Occupancy</h3><p>Phospho-S/T occupancy: the share of a protein’s copies carrying the phosphate, measured in 33 cancer cell lines (Koenig et al. 2026), median over the lines, unstimulated. Acetyl-K stoichiometry: the same for acetylation, measured in HeLa (Hansen et al. 2019).</p></div>' +
      '<div class="panel"><h3 id="def-tests">Tests</h3><dl>' +
      '<dt>Within protein</dt><dd>Mantel–Haenszel common odds ratio, each protein a stratum (its sites and its unmodified residues of the same type), with the Robins–Breslow–Greenland 95% CI and the continuity-corrected CMH test (PTM_Pipeline <code>common.mantel_haenszel_ci</code>; in the browser, the Methylarginine Site Explorer’s <code>mhTest</code>). It asks whether, within a protein, the modified residues lie in a context more often than its other residues, which a pooled test cannot separate from which proteins are modified.</dd>' +
      '<dt>Pooled</dt><dd>Fisher’s exact test with the Haldane–Anscombe odds ratio over the same counts (shown on hover).</dd>' +
      '<dt>Kinds of protein</dt><dd>One-sided hypergeometric test, fold enrichment (k/n ÷ M/N) against all reviewed human proteins.</dd>' +
      '<dt>Distributions</dt><dd>Mann–Whitney U, normal approximation with tie and continuity corrections.</dd>' +
      '<dt>Multiple testing</dt><dd>Benjamini–Hochberg across the rows of each chart; filled marks have q ≤ ' + PADJ + '. Effects within ' + C.minEffect() + '× either way are drawn grey (about the same) whatever their q.</dd>' +
      '<dt>Your selection</dt><dd>The same tests run again in the browser on the selected proteins and evidence, from the same per-protein counts; sites the selection does not keep are left out of both sides. The tick on each bar is the whole dataset. A selection’s q values are exploratory: they are not adjusted for the other selections tried.</dd></dl></div>' +
      '<div class="panel"><h3>How to read the figures</h3><p>Tiles and charts are the Methylarginine Site Explorer’s, drawn by its own code: a dot (or bar) at the odds ratio on a log scale, filled red where more often, slate where less, grey where about the same, hollow where not significant; a line through a dot is its 95% CI. Charts export as SVG, PNG and TSV.</p></div></div>';
  }

  // ---------------------------------------------------------------- dtIDR windows
  const inProtSel = () => (R.scoped && R.sel && R.sel.prots ? R.sel.prots : null);
  const selBanner = () => (inProtSel() ? '<p class="sel-banner">Showing the proteins of your selection: ' + esc(selectionText()) + ' · <a href="#' + parseHash().path + '">show all</a></p>' : '');
  const KCOLS = (D) => OV.kinds.map(([k, l]) => ({ key: k, label: l, num: true, get: (i) => D[k][i] }));
  async function flanksView(params) {
    view().innerHTML = '<p class="loading">Loading 121k windows…</p>';
    await Promise.all([loadProteins(), loadFlanks(), computeSelection()]);
    const n = FL.p.length;
    const st = { q: params.get('q') || '', sort: params.get('sort') || 'za', desc: params.get('asc') !== '1', page: +(params.get('page') || 1), chem: params.get('chem') ?? '', side: params.get('side') || '',
      pole: params.get('pole') || '', tf: params.get('tf') === '1', rbd: params.get('rbd') === '1', chap: params.get('chap') === '1', hist: params.get('hist') === '1', mr: params.get('mr') === '1', ck2: params.get('ck2') === '1' };
    view().innerHTML = '<h2>dtIDR windows</h2><p class="lede">Every window beside a folded domain, with both scores and their percentiles; the top tenth of each is its pole. Modifications are counted inside the window and never enter a score.' + docLink('unit') + '</p>' + selBanner() +
      '<div class="panel"><form class="filters" id="ff"><label>Gene or accession<input type="text" name="q" value="' + esc(st.q) + '" placeholder="e.g. NCL"></label>' +
      '<label>Sort by<select name="sort">' + [['za', 'acidic score'], ['zb', 'basic score'], ['both', 'both poles (lower of the two)'], ['lar', 'longest D/E run'], ['gene', 'gene']].map(([v, l]) => '<option value="' + v + '"' + (st.sort === v ? ' selected' : '') + '>' + l + '</option>').join('') + '</select></label>' +
      '<label>Pole<select name="pole"><option value="">any</option>' + [['a', 'acidic pole'], ['b', 'basic pole'], ['ab', 'both (bipolar)']].map(([v, l]) => '<option value="' + v + '"' + (st.pole === v ? ' selected' : '') + '>' + l + '</option>').join('') + '</select></label>' +
      '<label>Chemistry<select name="chem"><option value="">any</option>' + CHEM_SHORT.map((l, k) => '<option value="' + k + '"' + (st.chem === String(k) ? ' selected' : '') + '>' + l + '</option>').join('') + '</select></label>' +
      '<label>Side<select name="side"><option value="">either</option><option value="N"' + (st.side === 'N' ? ' selected' : '') + '>N-side</option><option value="C"' + (st.side === 'C' ? ' selected' : '') + '>C-side</option></select></label>' +
      [['tf', 'TF flank'], ['rbd', 'RNA-binding flank'], ['chap', 'histone chaperone or HMG'], ['mr', 'methyl-R inside'], ['ck2', 'CK2 site inside'], ['hist', 'include histone-like']].map(([k, l]) => '<label class="check"><input type="checkbox" name="' + k + '"' + (st[k] ? ' checked' : '') + '>' + l + '</label>').join('') +
      '</form></div><div class="plot-row"><div id="ft"></div><div class="panel" id="fplot"></div></div>';
    const q = st.q.trim().toUpperCase(), ps = inProtSel();
    const qAccs = q ? new Set(PR.acc.map((a, i) => ((PR.gene[i] || '').toUpperCase().startsWith(q) || a.startsWith(q) ? i : -1)).filter((i) => i >= 0)) : null;
    const rows = [];
    for (let i = 0; i < n; i++) {
      const f = FL.fl[i];
      if (!st.hist && (f & 4)) continue;
      if (ps && !ps.has(FL.p[i])) continue;
      if (qAccs && !qAccs.has(FL.p[i])) continue;
      if (st.chem !== '' && FL.chem[i] !== +st.chem) continue;
      if (st.side && FL.side[i] !== (st.side === 'N' ? 0 : 1)) continue;
      if (st.pole === 'a' && !(FL.pole[i] & 1)) continue;
      if (st.pole === 'b' && !(FL.pole[i] & 2)) continue;
      if (st.pole === 'ab' && FL.pole[i] !== 3) continue;
      if (st.tf && !(f & 1)) continue;
      if (st.rbd && !(f & 2)) continue;
      if (st.chap && !(f & 24)) continue;
      if (st.mr && !FL.mR[i]) continue;
      if (st.ck2 && !FL.ck2[i]) continue;
      rows.push(i);
    }
    const columns = [
      { key: 'gene', label: 'Gene', get: (i) => geneOf(FL.p[i]), html: (i) => protLink(FL.p[i], 'f=' + FL.fs[i] + '-' + FL.fe[i]) },
      { key: 'acc', label: 'Accession', get: (i) => PR.acc[FL.p[i]] },
      { key: 'side', label: 'Side', get: (i) => (FL.side[i] ? 'C' : 'N') },
      { key: 'fs', label: 'Window', get: (i) => FL.fs[i] + '-' + FL.fe[i], sort: (i) => FL.fs[i] },
      { key: 'dn', label: 'Adjacent domain', get: (i) => FL.domain_names[FL.dn[i]], html: (i) => esc((FL.domain_names[FL.dn[i]] || '').split(';')[0].slice(0, 48)) },
      { key: 'zb', label: 'Basic z', num: true, get: (i) => z(FL.zb[i]), html: (i) => zChip(FL.zb[i], 'basic') + (FL.pole[i] & 2 ? ' <span class="pole b" title="basic pole: top tenth">B</span>' : '') },
      { key: 'za', label: 'Acidic z', num: true, get: (i) => z(FL.za[i]), html: (i) => zChip(FL.za[i], 'acidic') + (FL.pole[i] & 1 ? ' <span class="pole a" title="acidic pole: top tenth">A</span>' : '') },
      { key: 'both', label: 'Lower of the two', num: true, get: (i) => (FL.zb[i] == null || FL.za[i] == null ? null : Math.min(FL.zb[i], FL.za[i]) / 100), html: (i) => fz(Math.min(FL.zb[i], FL.za[i])) },
      { key: 'chem', label: 'Chemistry', get: (i) => CHEM_SHORT[FL.chem[i]] || '' },
      { key: 'lar', label: 'Longest D/E run', num: true, get: (i) => FL.lar[i] },
    ].concat(KCOLS(FL), [{ key: 'flags', label: 'Class', get: (i) => flagText(FL.fl[i]) }]);
    const save = (over) => {
      Object.assign(st, over);
      const p = selParams();
      if (st.q) p.set('q', st.q);
      if (st.sort !== 'za') p.set('sort', st.sort);
      if (!st.desc) p.set('asc', '1');
      if (st.page > 1) p.set('page', st.page);
      for (const k of ['chem', 'side', 'pole']) if (st[k] !== '') p.set(k, st[k]);
      for (const k of ['tf', 'rbd', 'chap', 'mr', 'ck2', 'hist']) if (st[k]) p.set(k, '1');
      setParams('flanks', p);
    };
    const draw = () => table($('#ft'), { caption: 'dtIDR windows', columns, rows, sortKey: st.sort, desc: st.desc, page: st.page, csvName: 'dtidr_windows.tsv', onState: (o) => { save(o); draw(); $('#ft').scrollIntoView({ block: 'nearest' }); } });
    draw();
    const xs = FL.zb.map(z), ys = FL.za.map(z);
    let xr = [Infinity, -Infinity], yr = [Infinity, -Infinity];
    for (let i = 0; i < n; i++) { if (xs[i] != null) { xr[0] = Math.min(xr[0], xs[i]); xr[1] = Math.max(xr[1], xs[i]); } if (ys[i] != null) { yr[0] = Math.min(yr[0], ys[i]); yr[1] = Math.max(yr[1], ys[i]); } }
    xr = [Math.floor(xr[0] * 2) / 2, Math.ceil(xr[1] * 2) / 2]; yr = [Math.floor(yr[0] * 2) / 2, Math.ceil(yr[1] * 2) / 2];
    $('#fplot').innerHTML = '<h3 class="h-sm">The two poles</h3>' + density({ xs, ys, sel: rows.length < n ? rows : null, xLabel: 'basic dtIDR score (z)', yLabel: 'acidic dtIDR score (z)', xr, yr }) +
      '<p class="caveat">Binned density of all windows' + (rows.length < n ? '; the rows shown in orange' : '') + '. Most windows are 50 aa, so scores fall on a lattice.</p>';
    const form = $('#ff');
    form.addEventListener('change', (e) => { const t = e.target; if (t.type === 'checkbox') save({ [t.name]: t.checked, page: 1 }); else if (t.name === 'sort') save({ sort: t.value, desc: t.value !== 'gene', page: 1 }); else save({ [t.name]: t.value, page: 1 }); route(); });
    form.addEventListener('submit', (e) => e.preventDefault());
  }

  // ---------------------------------------------------------------- IDRs
  async function idrsView(params) {
    view().innerHTML = '<p class="loading">Loading 33k IDRs…</p>';
    await Promise.all([loadProteins(), loadIdrs(), computeSelection()]);
    const n = ID.p.length, b = (k) => params.get(k) === '1';
    const st = { q: params.get('q') || '', sort: params.get('sort') || 'dz', desc: params.get('asc') !== '1', page: +(params.get('page') || 1), minlen: +(params.get('minlen') || 0), term: params.get('term') || '',
      as: b('as'), rg: b('rg'), bp: b('bp'), mixed: b('mixed'), mr: b('mr'), glu: b('glu'), dt: b('dt'), ck2: b('ck2'), blocky: b('blocky'), sharp: b('sharp'), ncpr: params.get('ncpr') || '' };
    view().innerHTML = '<h2>IDRs</h2><p class="lede">Every MobiDB disordered region of ' + OV.settings.min_idr + ' or more residues: its charge and charge patterning (δ z), the change its observed modifications make (Δ NCPR, Δ δ), its sequence elements, its sites, and the best dtIDR window it overlaps.' + docLink('charge') + '</p>' + selBanner() +
      '<div class="panel"><form class="filters" id="if"><label>Gene or accession<input type="text" name="q" value="' + esc(st.q) + '" placeholder="e.g. NPM1"></label>' +
      '<label>Minimum length<input type="number" name="minlen" min="20" step="10" value="' + (st.minlen || '') + '" placeholder="20"></label>' +
      '<label>Net charge<select name="ncpr"><option value="">any</option><option value="neg"' + (st.ncpr === 'neg' ? ' selected' : '') + '>acidic (NCPR ≤ −0.1)</option><option value="pos"' + (st.ncpr === 'pos' ? ' selected' : '') + '>basic (NCPR ≥ +0.1)</option></select></label>' +
      '<label>Position<select name="term"><option value="">any</option>' + [['N', 'N-terminal'], ['C', 'C-terminal'], ['I', 'internal']].map(([v, l]) => '<option value="' + v + '"' + (st.term === v ? ' selected' : '') + '>' + l + '</option>').join('') + '</select></label>' +
      [['as', 'acidic stretch'], ['rg', 'RG/RGG cluster'], ['bp', 'basic patch'], ['mixed', 'acidic stretch and a basic patch or RG cluster'], ['blocky', 'blocky (δ z ≥ ' + OV.settings.blocky_z + ')'], ['sharp', 'phosphorylation sharpens its blocks more than random'], ['mr', 'methyl-R'], ['ck2', 'CK2 site'], ['glu', 'Glu site'], ['dt', 'overlaps a dtIDR window']]
        .map(([k, l]) => '<label class="check"><input type="checkbox" name="' + k + '"' + (st[k] ? ' checked' : '') + '>' + l + '</label>').join('') +
      '</form></div><div class="plot-row"><div id="it"></div><div class="panel" id="iplot"></div></div>';
    const q = st.q.trim().toUpperCase(), ps = inProtSel();
    const qAccs = q ? new Set(PR.acc.map((a, i) => ((PR.gene[i] || '').toUpperCase().startsWith(q) || a.startsWith(q) ? i : -1)).filter((i) => i >= 0)) : null;
    const rows = [];
    for (let i = 0; i < n; i++) {
      const len = ID.e[i] - ID.s[i] + 1;
      if (st.minlen && len < st.minlen) continue;
      if (ps && !ps.has(ID.p[i])) continue;
      if (qAccs && !qAccs.has(ID.p[i])) continue;
      if (st.term && ID.term[i] !== st.term) continue;
      if (st.ncpr === 'neg' && ID.ncpr[i] > -100) continue;
      if (st.ncpr === 'pos' && ID.ncpr[i] < 100) continue;
      if (st.as && !ID.nas[i]) continue;
      if (st.rg && !ID.nrg[i]) continue;
      if (st.bp && !ID.nbp[i]) continue;
      if (st.mixed && !(ID.nas[i] && (ID.nbp[i] || ID.nrg[i]))) continue;
      if (st.mr && !ID.mR[i]) continue;
      if (st.glu && !ID.glu[i]) continue;
      if (st.ck2 && !ID.ck2[i]) continue;
      if (st.blocky && !(ID.dz[i] >= OV.settings.blocky_z * 100)) continue;
      if (st.sharp && !(ID.dd_ph[i] != null && ID.rd_ph[i] != null && ID.dd_ph[i] > ID.rd_ph[i])) continue;
      if (st.dt && ID.zb[i] == null) continue;
      rows.push(i);
    }
    const f3 = (v) => (v == null ? '' : (v / 1000).toFixed(3));
    const fd = (v) => (v == null ? '' : (v >= 0 ? '+' : '−') + Math.abs(v / 1e6).toFixed(4));
    const columns = [
      { key: 'gene', label: 'Gene', get: (i) => geneOf(ID.p[i]), html: (i) => protLink(ID.p[i], 'f=' + ID.s[i] + '-' + ID.e[i]) },
      { key: 'acc', label: 'Accession', get: (i) => PR.acc[ID.p[i]] },
      { key: 'reg', label: 'IDR', get: (i) => ID.s[i] + '-' + ID.e[i], sort: (i) => ID.s[i] },
      { key: 'len', label: 'Length', num: true, get: (i) => ID.e[i] - ID.s[i] + 1 },
      { key: 'term', label: 'Position', get: (i) => ({ N: 'N-term', C: 'C-term', I: 'internal', W: 'whole' }[ID.term[i]]) },
      { key: 'ncpr', label: 'NCPR', num: true, get: (i) => ID.ncpr[i] / 1000, html: (i) => f3(ID.ncpr[i]) },
      { key: 'fcr', label: 'FCR', num: true, get: (i) => ID.fcr[i] / 1000, html: (i) => f3(ID.fcr[i]) },
      { key: 'dz', label: 'δ z (blocks)', num: true, get: (i) => z(ID.dz[i]), html: (i) => fz(ID.dz[i]) },
      { key: 'sz', label: 'SCD z', num: true, get: (i) => z(ID.sz[i]), html: (i) => fz(ID.sz[i]) },
      { key: 'dn', label: 'Δ NCPR (modified)', num: true, tip: 'Net charge per residue with every observed site modified, minus unmodified (idrkit.ptm charges)', get: (i) => (ID.nmod[i] ? ID.dn[i] / 1000 : null), html: (i) => (ID.nmod[i] ? f3(ID.dn[i]) : '') },
      { key: 'dd', label: 'Δ δ (modified)', num: true, tip: 'δ with every observed site modified, minus unmodified: above 0, the modifications sharpen the charge blocks', get: (i) => (ID.nmod[i] ? ID.dd[i] / 1e6 : null), html: (i) => (ID.nmod[i] ? fd(ID.dd[i]) : '') },
      { key: 'ddph', label: 'Δ δ phospho vs random', num: true, tip: 'δ change of the observed phosphosites minus the mean change of as many random S/T in the same IDR', get: (i) => (ID.dd_ph[i] != null && ID.rd_ph[i] != null ? (ID.dd_ph[i] - ID.rd_ph[i]) / 1e6 : null), html: (i) => (ID.dd_ph[i] != null && ID.rd_ph[i] != null ? fd(ID.dd_ph[i] - ID.rd_ph[i]) : '') },
      { key: 'lde', label: 'Longest D/E run', num: true, get: (i) => ID.lde[i] },
      { key: 'nas', label: 'Acidic stretches', num: true, get: (i) => ID.nas[i] },
      { key: 'nbp', label: 'Basic patches', num: true, get: (i) => ID.nbp[i] },
      { key: 'nrg', label: 'RG clusters', num: true, get: (i) => ID.nrg[i] },
    ].concat(KCOLS(ID), [
      { key: 'zb', label: 'dtIDR basic z', num: true, get: (i) => z(ID.zb[i]), html: (i) => zChip(ID.zb[i], 'basic') },
      { key: 'za', label: 'dtIDR acidic z', num: true, get: (i) => z(ID.za[i]), html: (i) => zChip(ID.za[i], 'acidic') }]);
    const save = (over) => {
      Object.assign(st, over);
      const p = selParams();
      if (st.q) p.set('q', st.q);
      if (st.sort !== 'dz') p.set('sort', st.sort);
      if (!st.desc) p.set('asc', '1');
      if (st.page > 1) p.set('page', st.page);
      if (st.minlen) p.set('minlen', st.minlen);
      if (st.term) p.set('term', st.term);
      if (st.ncpr) p.set('ncpr', st.ncpr);
      for (const k of ['as', 'rg', 'bp', 'mixed', 'blocky', 'sharp', 'mr', 'ck2', 'glu', 'dt']) if (st[k]) p.set(k, '1');
      setParams('idrs', p);
    };
    const draw = () => table($('#it'), { caption: 'IDRs', columns, rows, sortKey: st.sort, desc: st.desc, page: st.page, csvName: 'idrs.tsv', onState: (o) => { save(o); draw(); $('#it').scrollIntoView({ block: 'nearest' }); } });
    draw();
    $('#iplot').innerHTML = '<h3 class="h-sm">Charge of every IDR</h3>' + density({ xs: ID.fcr.map((v) => v / 1000), ys: ID.ncpr.map((v) => v / 1000), sel: rows.length < n ? rows : null, xLabel: 'fraction of charged residues', yLabel: 'net charge per residue', xr: [0, 1], yr: [-1, 1], bins: 50 }) +
      '<p class="caveat">K and R +1, D and E −1, His 0. Rows shown in orange.</p>';
    $('#if').addEventListener('change', (e) => { const t = e.target; if (t.type === 'checkbox') save({ [t.name]: t.checked, page: 1 }); else if (t.name === 'minlen') save({ minlen: +t.value || 0, page: 1 }); else save({ [t.name]: t.value, page: 1 }); route(); });
    $('#if').addEventListener('submit', (e) => e.preventDefault());
  }

  // ---------------------------------------------------------------- protein page
  const EMBOSS = { K: 1, R: 1, H: 0.5, D: -1, E: -1 };
  function chargeProfile(seq, bins) {
    const n = seq.length, w = 2, q = Array.from(seq, (c) => EMBOSS[c] || 0), m = new Float32Array(n);
    for (let i = 0; i < n; i++) { let s = 0, k = 0; for (let j = Math.max(0, i - w); j <= Math.min(n - 1, i + w); j++) { s += q[j]; k++; } m[i] = s / k; }
    const out = [];
    for (let b = 0; b < bins; b++) { const a = Math.floor(b * n / bins), e = Math.max(a + 1, Math.floor((b + 1) * n / bins)); let s = 0; for (let i = a; i < e; i++) s += m[i]; out.push(s / (e - a)); }
    return out;
  }
  async function proteinView(key, params, token) {
    await loadProteins();
    const K = key.toUpperCase(), p = byAcc.has(K) ? byAcc.get(K) : byGene.get(K);
    if (p == null) { view().innerHTML = '<div class="panel">No reviewed human protein called <b>' + esc(key) + '</b>.</div>'; return; }
    const acc = PR.acc[p];
    const [sh] = await Promise.all([loadShard(acc), loadIdrs()]);
    if (token !== routeToken) return;
    const Rp = sh.prot[acc], names = sh.names, seq = Rp.seq, n = seq.length, S = Rp.s || {};
    const hl = (params.get('f') || '').match(/^(\d+)-(\d+)$/);
    const W = 1000, L = 150, PW = W - L - 10;
    const x = (pos) => L + (pos - 1) / Math.max(1, n - 1) * PW, wd = (s, e) => Math.max(1.5, (e - s + 1) / n * PW);
    let y = 8;
    const rows = [];
    const row = (label, h, inner) => { rows.push('<text class="lab" x="' + (L - 8) + '" y="' + (y + h / 2 + 4) + '" text-anchor="end">' + esc(label) + '</text>' + inner(y, h)); y += h + 8; };
    const seg = (s, e, yy, h, style, tip) => '<rect x="' + x(s).toFixed(1) + '" y="' + yy + '" width="' + wd(s, e).toFixed(1) + '" height="' + h + '" rx="2" style="' + style + '" data-s="' + s + '" data-e="' + e + '" data-tip="' + esc(tip) + '"/>';
    const name = (k) => (names[k] ? names[k][1] : '');
    row('Domains', 16, (yy, h) => (Rp.dom || []).map((d) => seg(d[0], d[1], yy, h, 'fill:var(--domain);stroke:var(--domain-ink);stroke-width:0.6', '<b>' + esc(name(d[2])) + '</b> ' + d[0] + '–' + d[1])).join(''));
    if ((Rp.fam || []).length) row('Families', 8, (yy, h) => Rp.fam.map((d) => seg(d[0], d[1], yy, h, 'fill:var(--family)', esc(name(d[2])) + ' ' + d[0] + '–' + d[1])).join(''));
    row('Disorder (MobiDB)', 8, (yy, h) => '<line x1="' + L + '" x2="' + (L + PW) + '" y1="' + (yy + h / 2) + '" y2="' + (yy + h / 2) + '" style="stroke:var(--line-strong)"/>' + (Rp.dis || []).map((d) => seg(d[0], d[1], yy, h, 'fill:var(--dis)', 'Disordered ' + d[0] + '–' + d[1])).join('') +
      (Rp.blocky || []).map((d) => seg(d[0], d[1], yy - 2, h + 4, 'fill:none;stroke:var(--ink);stroke-width:1.2', '<b>Blocky IDR</b> ' + d[0] + '–' + d[1] + ' (δ z ≥ ' + OV.settings.blocky_z + ')')).join(''));
    const flankRow = (label, idx, pctIdx, pole, bit) => row(label, 14, (yy, h) => (Rp.fl || []).map((f) => {
      const top = (f[10] || 0) & bit, pc = f[pctIdx], a = top ? 0.9 : pc == null ? 0.05 : Math.max(0.05, (pc / 1000 - 0.5) * 1.2);
      return seg(f[1], f[2], yy, h, 'fill:var(--' + pole + ');fill-opacity:' + a.toFixed(2) + ';stroke:var(--' + pole + ');stroke-opacity:' + (top ? 1 : 0.5) + ';stroke-width:' + (top ? 1.2 : 0.5),
        '<b>' + (f[0] ? 'C' : 'N') + '-side window ' + f[1] + '–' + f[2] + '</b> (domain ' + f[3] + '–' + f[4] + ')<br>' + pole + ' z ' + fz(f[idx]) + ', percentile ' + (pc / 10).toFixed(1) + (top ? ' · <b>' + pole + ' pole</b>' : '') + (f[7] >= 0 ? '<br>' + CHEM_SHORT[f[7]] : ''));
    }).join(''));
    flankRow('dtIDR, basic pole', 5, 8, 'basic', 2);
    flankRow('dtIDR, acidic pole', 6, 9, 'acidic', 1);
    const bins = Math.min(n, 400), prof = chargeProfile(seq, bins);
    row('Charge (5-aa mean)', 40, (yy, h) => {
      const mid = yy + h / 2, bw = PW / bins;
      let s = '<line x1="' + L + '" x2="' + (L + PW) + '" y1="' + mid + '" y2="' + mid + '" style="stroke:var(--line-strong)"/>';
      prof.forEach((v, b) => { if (!v) return; const hh = Math.abs(v) * h / 2; s += '<rect x="' + (L + b * bw).toFixed(1) + '" y="' + (v > 0 ? mid - hh : mid).toFixed(1) + '" width="' + (bw + 0.2).toFixed(2) + '" height="' + hh.toFixed(1) + '" style="fill:var(--' + (v > 0 ? 'basic' : 'acidic') + ')"/>'; });
      return s;
    });
    row('Acidic stretches', 8, (yy, h) => (Rp.as || []).map((d) => seg(d[0], d[1], yy, h, 'fill:var(--acidic)', 'Acidic stretch ' + d[0] + '–' + d[1])).join(''));
    row('Basic patches', 8, (yy, h) => (Rp.bp || []).map((d) => seg(d[0], d[1], yy, h, 'fill:var(--basic)', 'Basic patch ' + d[0] + '–' + d[1] + ' (net ≥ +3 in 10 aa)')).join(''));
    row('RG/RGG clusters', 8, (yy, h) => (Rp.rg || []).map((d) => seg(d[0], d[1], yy, h, 'fill:var(--rg);stroke:var(--rg-ink);stroke-width:0.6', 'RG/RGG cluster ' + d[0] + '–' + d[1])).join(''));
    const KIN = OV.kinase_rows;
    const kinText = (v) => (v.length > 4 && v[3] != null ? '<br>CK2 percentile ' + v[3] + '; motif matches: ' + (KIN.filter((_, b) => v[4] & (1 << b)).map((r) => r[1]).join(', ') || 'none ≥ ' + OV.kinase_match) + '<br>best-scoring kinase ' + esc(B.top_kinases[v[5]]) + ' (' + v[6] + ')' + (v[7] != null ? '<br>occupancy ' + (v[7] / 10).toFixed(1) + '% (median of cell lines)' : '') : '');
    const ticks = (label, list, color, tip, tall = () => false, hollow = null) => row(label, 14, (yy, h) => (list || []).map((v) => {
      const pos = v[0], xx = x(pos).toFixed(1), t = tall(v);
      return '<line x1="' + xx + '" x2="' + xx + '" y1="' + (t ? yy - 2 : yy + 3) + '" y2="' + (yy + h) + '" style="stroke:' + color(v) + ';stroke-width:' + (t ? 2.2 : 1.4) + '" data-s="' + pos + '" data-e="' + pos + '" data-tip="' + esc(tip(v)) + '"/>' +
        (hollow ? '<circle cx="' + xx + '" cy="' + (yy - 1) + '" r="3" style="fill:' + (hollow(v) ? 'var(--surface)' : color(v)) + ';stroke:' + color(v) + '" data-tip="' + esc(tip(v)) + '"/>' : '');
    }).join(''));
    const src = (v) => srcNames(v[1]) + (v[2] ? ', PTMeXchange ' + CONF[v[2]] : '');
    const tip = (what) => (v) => '<b>' + seq[v[0] - 1] + v[0] + ' ' + what + '</b><br>' + src(v);
    // rows by residue: arginine marks, lysine marks, serine/threonine/tyrosine marks, glutamate (the order the palette was checked in)
    ticks('Methyl-R', S.mR, () => 'var(--meth)', (v) => tip('methyl-R')(v) + (v[3] ? '<br>' + esc(v[3]) : ''));
    if ((S.cit || []).length) ticks('Citrulline', S.cit, () => 'var(--cit)', tip('citrulline'));
    ticks('Acetyl-K', S.kac, () => 'var(--kac)', (v) => tip('acetyl')(v) + (v[3] != null && typeof v[3] === 'number' ? '<br>stoichiometry ' + (100 * v[3]).toFixed(3) + '% (Hansen 2019)' : v[3] ? '<br>' + esc(v[3]) : ''), (v) => v[2] === 3);
    if ((S.kacyl || []).length) ticks('Acidic acyl-K', S.kacyl, () => 'var(--kac)', (v) => tip('acidic acyl')(v) + (v[3] ? '<br>' + esc(v[3]) : ''));
    if ((S.ub || []).length) ticks('Ubiquitin-K', S.ub, () => 'var(--phospho)', tip('ubiquitin'), (v) => v[2] === 3);
    if ((S.sumo || []).length) ticks('SUMO-K', S.sumo, () => 'var(--sumo)', tip('SUMO'), (v) => v[2] === 3);
    if ((S.kme || []).length) ticks('Methyl-K', S.kme, () => 'var(--kme)', tip('methyl-K'));
    ticks('Phospho-S/T (CK2 tall)', S.ph, (v) => ((S.ck2 || []).some((c) => c[0] === v[0]) ? 'var(--ck2)' : 'var(--phospho)'), (v) => tip('phospho')(v) + kinText(v), (v) => (S.ck2 || []).some((c) => c[0] === v[0]));
    if ((S.ck2m || []).length) ticks('Strong CK2 motifs, unobserved', S.ck2m, () => 'var(--ck2)', (v) => '<b>' + seq[v[0] - 1] + v[0] + '</b>: CK2 percentile ' + v[1] + ', no phosphosite reported', () => false, () => true);
    if ((S.py || []).length) ticks('Phospho-Y', S.py, () => 'var(--phospho)', tip('phospho-Y'));
    if ((S.adpr || []).length) ticks('ADP-ribosyl', S.adpr, () => 'var(--adpr)', tip('ADP-ribosyl'));
    if ((S.glu || []).length) ticks('Glutamylation', S.glu, () => 'var(--glu)', (v) => tip('glutamylation')(v) + (v[3] ? '<br>' + esc(v[3]) : ''));
    const axisY = y;
    let axis = '<line x1="' + L + '" x2="' + (L + PW) + '" y1="' + axisY + '" y2="' + axisY + '" style="stroke:var(--muted)"/>';
    const step = [10, 20, 50, 100, 200, 250, 500, 1000, 2000, 5000].find((s) => n / s <= 10) || 10000;
    for (let t = 0; t <= n; t += step) { const xx = x(Math.max(1, t)); axis += '<line x1="' + xx + '" x2="' + xx + '" y1="' + axisY + '" y2="' + (axisY + 4) + '" style="stroke:var(--muted)"/><text x="' + xx + '" y="' + (axisY + 16) + '" text-anchor="middle" font-size="11">' + Math.max(1, t) + '</text>'; }
    const H = axisY + 24;
    const band = hl ? '<rect x="' + x(+hl[1]).toFixed(1) + '" y="0" width="' + wd(+hl[1], +hl[2]).toFixed(1) + '" height="' + axisY + '" style="fill:var(--accent);fill-opacity:0.12;stroke:var(--accent);stroke-opacity:0.6"/>' : '';
    const svg = '<svg class="map" id="pmap" viewBox="0 0 ' + W + ' ' + H + '" role="img" aria-label="Protein map of ' + esc(geneOf(p)) + '">' + band + rows.join('') + axis + '</svg>';
    const myIdrs = [];
    for (let i = 0; i < ID.p.length; i++) if (ID.p[i] === p) myIdrs.push(i);
    const legend = [['var(--domain)', 'InterPro domain'], ['var(--dis)', 'disordered (outlined: blocky)'], ['var(--basic)', 'basic: pole (solid: top tenth), charge, patch'], ['var(--acidic)', 'acidic: pole (solid: top tenth), charge, stretch'],
      ['var(--rg)', 'RG/RGG cluster'], ['var(--meth)', 'methyl-R'], ['var(--cit)', 'citrulline'], ['var(--kac)', 'acetyl / acyl-K'], ['var(--phospho)', 'ubiquitin-K, phospho'], ['var(--sumo)', 'SUMO-K'], ['var(--kme)', 'methyl-K'],
      ['var(--ck2)', 'CK2 phospho (hollow: strong motif, unobserved)'], ['var(--adpr)', 'ADP-ribosyl'], ['var(--glu)', 'glutamylation']].map(([c, l]) => '<span><i style="background:' + c + '"></i>' + l + '</span>').join('');
    view().innerHTML = '<div class="prot-head"><h2>' + esc(geneOf(p)) + ' <span class="muted acc">' + esc(acc) + '</span></h2><p class="pmeta">' + esc(PR.name[p]) + ' · ' + fmt(n) + ' aa · ' + PR.dis[p] + '% disordered</p>' +
      '<div class="links"><a href="' + PTM_SITE + '#' + acc + '">Methylarginine Site Explorer</a><a href="https://www.uniprot.org/uniprotkb/' + acc + '/entry">UniProt</a><a href="https://alphafold.ebi.ac.uk/entry/' + acc + '">AlphaFold</a><a href="https://mobidb.org/' + acc + '">MobiDB</a></div></div>' +
      '<div class="panel"><div class="map-scroll">' + svg + '</div><div class="legend">' + legend + '</div><p class="caveat">Hover a mark for its sources and scores; click a segment or site to read its sequence. Tall marks: Gold PTMeXchange sites, or CK2-attributed phosphosites.</p><div id="segment"></div></div>' +
      '<div id="pfl"></div><div id="pidr" class="mt"></div>';
    const showSeg = (s, e, label) => {
      s = Math.max(1, s - (s === e ? 15 : 0)); e = Math.min(n, e + (s === e ? 15 : 0));
      const marks = new Set();
      for (const k of Object.keys(S)) if (k !== 'ck2m') for (const v of S[k]) marks.add(v[0]);
      let html = '';
      for (let i = s; i <= e; i++) { const c = seq[i - 1], cls = ('KR'.includes(c) ? 'b' : 'DE'.includes(c) ? 'a' : '') + (marks.has(i) ? ' s' : ''); html += cls ? '<span class="' + cls.trim() + '" title="' + c + i + '">' + c + '</span>' : c; if ((i - s + 1) % 10 === 0) html += ' '; }
      $('#segment').innerHTML = '<h3 class="h-sm">' + esc(label || '') + ' ' + s + '–' + e + '</h3><div class="seq">' + html + '</div><p class="caveat">Basic residues in blue, acidic in red; modified residues underlined.</p>';
    };
    $('#pmap').addEventListener('click', (e) => { const t = e.target.closest('[data-s]'); if (t) showSeg(+t.dataset.s, +t.dataset.e, ''); });
    if (hl) showSeg(+hl[1], +hl[2], 'Selected');
    table($('#pfl'), { caption: 'dtIDR windows of ' + geneOf(p), rows: (Rp.fl || []).map((f, k) => k), sortKey: 'fs', desc: false, page: 1, columns: [
      { key: 'fs', label: 'Window', get: (k) => Rp.fl[k][1] + '-' + Rp.fl[k][2], sort: (k) => Rp.fl[k][1], html: (k) => '<a href="#' + acc + '?f=' + Rp.fl[k][1] + '-' + Rp.fl[k][2] + '">' + Rp.fl[k][1] + '–' + Rp.fl[k][2] + '</a>' },
      { key: 'side', label: 'Side', get: (k) => (Rp.fl[k][0] ? 'C' : 'N') }, { key: 'dom', label: 'Domain', get: (k) => Rp.fl[k][3] + '-' + Rp.fl[k][4] },
      { key: 'zb', label: 'Basic z', num: true, get: (k) => z(Rp.fl[k][5]), html: (k) => zChip(Rp.fl[k][5], 'basic') + ((Rp.fl[k][10] || 0) & 2 ? ' <span class="pole b">B</span>' : '') },
      { key: 'pb', label: 'Basic percentile', num: true, get: (k) => Rp.fl[k][8] / 10 },
      { key: 'za', label: 'Acidic z', num: true, get: (k) => z(Rp.fl[k][6]), html: (k) => zChip(Rp.fl[k][6], 'acidic') + ((Rp.fl[k][10] || 0) & 1 ? ' <span class="pole a">A</span>' : '') },
      { key: 'pa', label: 'Acidic percentile', num: true, get: (k) => Rp.fl[k][9] / 10 }, { key: 'chem', label: 'Chemistry', get: (k) => CHEM_SHORT[Rp.fl[k][7]] || '' }] });
    const f3 = (v) => (v / 1000).toFixed(3);
    table($('#pidr'), { caption: 'IDRs of ' + geneOf(p) + ' (≥' + OV.settings.min_idr + ' aa)', rows: myIdrs, sortKey: 's', desc: false, page: 1, columns: [
      { key: 's', label: 'IDR', get: (i) => ID.s[i] + '-' + ID.e[i], sort: (i) => ID.s[i], html: (i) => '<a href="#' + acc + '?f=' + ID.s[i] + '-' + ID.e[i] + '">' + ID.s[i] + '–' + ID.e[i] + '</a>' },
      { key: 'len', label: 'Length', num: true, get: (i) => ID.e[i] - ID.s[i] + 1 }, { key: 'ncpr', label: 'NCPR', num: true, get: (i) => f3(ID.ncpr[i]) },
      { key: 'dz', label: 'δ z (blocks)', num: true, get: (i) => z(ID.dz[i]), html: (i) => fz(ID.dz[i]) },
      { key: 'dn', label: 'Δ NCPR (modified)', num: true, get: (i) => (ID.nmod[i] ? f3(ID.dn[i]) : '') },
      { key: 'nas', label: 'Acidic stretches', num: true, get: (i) => ID.nas[i] }, { key: 'nbp', label: 'Basic patches', num: true, get: (i) => ID.nbp[i] },
      { key: 'mR', label: 'methyl-R', num: true, get: (i) => ID.mR[i] }, { key: 'ph', label: 'phospho', num: true, get: (i) => ID.ph[i] }, { key: 'cit', label: 'citrulline', num: true, get: (i) => ID.cit[i] }] });
    document.title = geneOf(p) + ' · IDR Explorer';
  }

  // ---------------------------------------------------------------- sources
  async function aboutView() {
    const c = B.counts;
    const inputs = Object.entries(B.inputs).map(([k, v]) => '<tr><td>' + esc(k) + '</td><td>' + esc(v.repo) + '</td><td><code>' + esc(v.path) + '</code></td><td><code>' + esc(v.md5.slice(0, 10)) + '…</code></td></tr>').join('');
    const commit = (repo, sha, pub) => (sha ? (pub ? '<a href="https://github.com/Shechterlab/' + repo + '/commit/' + sha + '"><code>' + sha.slice(0, 9) + '</code></a>' : '<code>' + sha.slice(0, 9) + '</code>') : 'not recorded');
    const kinds = c.by_kind;
    view().innerHTML = '<div class="about"><h2>Sources</h2><div class="panel"><p>Generated, not edited: <code>Analysis/scripts/build_idr_explorer.py</code> in MicroDomains, with PTM_Pipeline, glu_proteomics and this repository cloned beside it. Definitions live in one place each: sequence rules in <a href="' + GH.site + 'idrkit">idrkit</a> (this repository), regions and sites in MicroDomains’ <code>idr_universe.py</code>, staged public layers in PTM_Pipeline. The Overview’s charts and tests are the Methylarginine Site Explorer’s own code (<a href="' + GH.site + 'tools/sync-ptm-core.mjs">tools/sync-ptm-core.mjs</a>).</p>' +
      '<dl><dt>Built</dt><dd>' + esc(B.generated) + ' from PTM_Pipeline ' + commit('PTM_Pipeline', B.commits.PTM_Pipeline, true) + ', idr-explorer ' + commit('idr-explorer', B.commits['idr-explorer'], true) + ', MicroDomains ' + commit('', B.commits.MicroDomains) + ' and glu_proteomics ' + commit('', B.commits.glu_proteomics) + ' (private until publication).</dd>' +
      '<dt>Reference, disorder, domains</dt><dd>UniProt 2026_03 reviewed human (CC BY 4.0); MobiDB (Piovesan et al., Nucleic Acids Res 53:D495, 2025; CC BY 4.0); InterPro 108.0 (CC0).</dd>' +
      '<dt>Methylarginines</dt><dd>PTM_Pipeline release ' + esc(OV.settings.release) + ' (' + fmt(kinds.mR) + ' with PTMeXchange’s): DeAngelo, Silverstein & Shechter, Cell Chem Biol 2026, doi:10.1016/j.chembiol.2026.09.011 (CC BY-NC-SA 4.0).</dd>' +
      '<dt>Phosphorylation</dt><dd>' + fmt(kinds.ph) + ' phospho-S/T and ' + fmt(kinds.py) + ' phospho-Y: iPTMnet (Huang et al., Nucleic Acids Res 46:D542, 2018; CC BY-NC-SA 4.0), UniProt 2026_03 MOD_RES with its kinases (' + fmt(kinds.ck2) + ' CK2-attributed), and the class I sites of the phosphoproteome atlas of 33 cell lines (Koenig et al., Nat Struct Mol Biol 33:1561, 2026, doi:10.1038/s41594-026-01877-6; CC BY 4.0) with its occupancy. No PhosphoSitePlus.</dd>' +
      '<dt>Lysine marks</dt><dd>Acetyl (' + fmt(kinds.kac) + '; stoichiometry from Hansen et al., Nat Commun 10:1055, 2019, doi:10.1038/s41467-019-09024-0; CC BY 4.0), acidic acyl (' + fmt(kinds.kacyl) + '), ubiquitin (' + fmt(kinds.ub) + '), SUMO (' + fmt(kinds.sumo) + ') and methyl (' + fmt(kinds.kme) + ') lysines from iPTMnet, UniProt and PTMeXchange’s reprocessed human builds (PRIDE PXD079218, PXD068989, PXD068489, PXD076372; CC0), graded Gold, Silver or Bronze.</dd>' +
      '<dt>Citrullination</dt><dd>' + fmt(kinds.cit) + ' sites: the citrullinome of Rebak et al., Nat Struct Mol Biol 31:977, 2024 (doi:10.1038/s41594-024-01214-9; CC BY 4.0) and UniProt.</dd>' +
      '<dt>ADP-ribosylation</dt><dd>' + fmt(kinds.adpr) + ' sites: Hendriks et al., Nat Commun 12:5893, 2021 (doi:10.1038/s41467-021-26172-4; CC BY 4.0) and UniProt.</dd>' +
      '<dt>Glutamylation</dt><dd>' + fmt(kinds.glu) + ' sites: the TTLL4 re-analysis of PRIDE PXD065251 (lab analysis), its published sites and UniProt’s polyglutamylation.</dd>' +
      '<dt>Kinase motifs</dt><dd>Johnson et al., Nature 613:759, 2023 (doi:10.1038/s41586-022-05575-3; CC BY 4.0), scored with kinase-library 1.8.0 (CC BY-NC-SA 3.0): ' + fmt(c.kinase_scored) + ' serines and threonines.</dd>' +
      '<dt>Partners</dt><dd>STRING v12.5 physical (CC BY 4.0), Complex Portal (CC0), IntAct binding regions (CC BY 4.0), staged in PTM_Pipeline.</dd>' +
      '<dt>Analyses</dt><dd>The dtIDR scores and the Overview’s statistics are the Shechter lab’s MicroDomains analyses (private until publication), built from the public inputs above.</dd></dl></div>' +
      '<div class="panel"><h3 class="h-sm">Inputs</h3><div class="table-wrapper" tabindex="0" role="region" aria-label="Inputs and checksums"><table class="data"><thead><tr><th>Input</th><th>Repository</th><th>Path</th><th>MD5</th></tr></thead><tbody>' + inputs + '</tbody></table></div></div>' +
      '<div class="panel"><h3 class="h-sm">Read before quoting</h3><ul><li>Modification sites are reported sites: they track study effort as well as biology. The within-protein tests remove differences between proteins, not between regions a method sees well.</li>' +
      '<li>Kinase motifs are predictions from sequence. Charge models are one modification per residue at pH 7.4.</li><li>Scores live on a lattice (most windows are 50 aa); ties are expected.</li><li>Glutamylation sites from MS favour peptides beside K/R; most are region-level.</li></ul></div></div>';
  }

  // ---------------------------------------------------------------- search and theme
  const input = $('#q'), sug = $('#suggest');
  const matches = (t) => {
    const T = t.trim().toUpperCase();
    if (!T || !PR) return [];
    const out = [];
    if (byAcc.has(T)) out.push(byAcc.get(T));
    const genes = [];
    for (let i = 0; i < PR.acc.length; i++) if ((PR.gene[i] || '').toUpperCase().startsWith(T) && !out.includes(i)) genes.push(i);
    genes.sort((a, b) => PR.gene[a].length - PR.gene[b].length || PR.gene[a].localeCompare(PR.gene[b]));
    out.push(...genes.slice(0, 8 - out.length));
    for (let i = 0; i < PR.acc.length && out.length < 8; i++) if (PR.acc[i].startsWith(T) && !out.includes(i)) out.push(i);
    return out;
  };
  input.addEventListener('focus', () => { loadProteins(); }, { once: true });
  input.addEventListener('input', async () => {
    await loadProteins();
    const m = matches(input.value);
    sug.hidden = !m.length;
    sug.innerHTML = m.map((i) => '<li><a href="#' + PR.acc[i] + '">' + esc(PR.gene[i] || PR.acc[i]) + ' <span class="muted">' + PR.acc[i] + ' · ' + esc(PR.name[i].slice(0, 40)) + '</span></a></li>').join('');
  });
  input.addEventListener('keydown', async (e) => {
    if (e.key === 'Enter') { await loadProteins(); const m = matches(input.value); if (m.length) location.hash = PR.acc[m[0]]; sug.hidden = true; }
    if (e.key === 'ArrowDown' && !sug.hidden) { e.preventDefault(); const a = sug.querySelector('a'); if (a) a.focus(); }
    if (e.key === 'Escape') sug.hidden = true;
  });
  sug.addEventListener('click', () => { sug.hidden = true; input.value = ''; });
  document.addEventListener('click', (e) => { if (!e.target.closest('.search')) sug.hidden = true; });
  // the explorer's colour tokens (--t-*, its TOKENS), which its stylesheet reads: set as it sets them, with the theme
  C.applyTokens();
  try { matchMedia('(prefers-color-scheme: dark)').addEventListener('change', () => C.applyTokens()); } catch (e) { /* old browsers */ }
  $('#theme-btn').addEventListener('click', () => {
    const cur = document.documentElement.getAttribute('data-theme') || (matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light');
    const next = cur === 'dark' ? 'light' : 'dark';
    document.documentElement.setAttribute('data-theme', next);
    try { localStorage.setItem('idr-theme', next); } catch (e) { /* private mode */ }
    C.applyTokens(next);
  });
  loadBase().then(() => { $('#built').textContent = 'Data built ' + B.generated + '.'; }).catch(() => {});
  route();
})();
