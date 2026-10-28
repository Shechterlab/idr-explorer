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
  let FAM = null;
  const loadFam = async () => (FAM = FAM || await getJSON('data/families.json'));
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
  const z = (v) => (v == null ? null : v / 100);
  const fz = (v) => (v == null ? '' : (v >= 0 ? '+' : '−') + Math.abs(v / 100).toFixed(2));
  const zChip = (v, pole) => {
    if (v == null) return '<span class="muted">–</span>';
    const a = Math.max(0, Math.min(1, (v / 100) / 3));
    return '<span class="z" style="background:color-mix(in srgb, var(--' + pole + ') ' + Math.round(a * 55) + '%, transparent)">' + fz(v) + '</span>';
  };
  const f3 = (v) => (v == null ? '' : (v / 1000).toFixed(3));
  // Das & Pappu's diagram of states (PNAS 110:13392, 2013), the classes localCIDER assigns (get_phasePlotRegion): the
  // field's names, a plain name, and the shape each class predicts for a disordered chain
  const REGION = ['', 'weak polyampholyte or polyelectrolyte', 'Janus sequence', 'strong polyampholyte', 'strong polyelectrolyte, negative', 'strong polyelectrolyte, positive'];
  const REGION_SHORT = ['', 'R1 weakly charged', 'R2 in between', 'R3 strongly charged, mixed', 'R4 strongly negative', 'R5 strongly positive'];
  const REGION_SHAPE = ['', 'compact (globule or tadpole)', 'compact or expanded, depending on context', 'coil or hairpin, depending on how + and − are mixed (κ)', 'expanded (swollen coil)', 'expanded (swollen coil)'];
  const REGION_RULE = ['', 'fewer than a quarter of residues charged (FCR < 0.25)', 'a quarter to a third charged (FCR 0.25–0.35)', 'over a third charged, balanced (FCR > 0.35, |net| < 0.35)',
    'over a third of residues negative (f− > 0.35)', 'over a third of residues positive (f+ > 0.35)'];
  const REGION_DEF = REGION_RULE.map((r, k) => (k ? r + '; predicted shape: ' + REGION_SHAPE[k] : ''));
  const dpTip = (r) => '<b>Das–Pappu ' + REGION_SHORT[r] + '</b><br>' + REGION[r] + ': ' + REGION_RULE[r] + '<br>Predicted shape: ' + REGION_SHAPE[r] + '<br><span class="muted">Das & Pappu 2013, as localCIDER assigns it</span>';
  const dpChip = (r, short) => (r ? '<span class="dp dp' + r + '" data-tip="' + esc(dpTip(r)) + '">' + (short ? 'R' + r : esc(REGION_SHORT[r])) + '</span>' : '<span class="muted">–</span>');
  const fsign = (v, d = 0) => (v == null || !Number.isFinite(v) ? '—' : (v > 0 ? '+' : v < 0 ? '−' : '') + Math.abs(v).toFixed(d));
  // Charge in words (idrkit.charge_type): net charge per residue first, then the share of residues charged.
  const CT = ['', 'Few charges', 'Mixed charge', 'Acidic', 'Basic'];
  const CT_DEF = ['', 'net charge between −0.10 and +0.10 per residue, and fewer than a quarter of residues charged', 'net charge between −0.10 and +0.10 per residue, with a quarter or more of residues charged (positive and negative balanced)',
    'net charge −0.10 per residue or below (5 or more negative per 50 residues)', 'net charge +0.10 per residue or above (5 or more positive per 50 residues)'];
  const ctOf = (fcr, ncpr) => (ncpr <= -0.1 ? 3 : ncpr >= 0.1 ? 4 : fcr >= 0.25 ? 2 : 1);
  const ctChip = (t) => (t ? '<span class="ct ct' + t + '" data-tip="' + esc('<b>' + CT[t] + '</b>: ' + CT_DEF[t]) + '">' + CT[t] + '</span>' : '<span class="muted">–</span>');
  // D/E (red) and K/R (blue) as bars over 0-50% of residues, and the net charge per residue
  const chargeCell = (fm, fp, net, t) => '<span class="cbar" data-tip="' + esc((t ? '<b>' + CT[t] + '</b><br>' : '') + 'D/E ' + Math.round(100 * fm) + '% · K/R ' + Math.round(100 * fp) + '% of residues<br>net charge ' + fsign(net, 2) + ' per residue') + '"><i class="m" style="width:' + Math.min(100, 200 * fm).toFixed(0) + '%"></i><i class="p" style="width:' + Math.min(100, 200 * fp).toFixed(0) + '%"></i></span><span class="cnet">' + fsign(net, 2) + '</span>';
  const dpChips = (bits) => [1, 2, 3, 4, 5].filter((r) => bits & (1 << (r - 1))).map((r) => dpChip(r, true)).join('') || '<span class="muted">–</span>';

  // ---------------------------------------------------------------- selection (the explorer's F and R)
  // Proteins: gene list, regions they have (all must hold), classes (any of those ticked). Evidence: the sources a site
  // must come from (any of them), and the lowest PTMeXchange grade for a site PTMeXchange alone reports.
  const HAS_LABEL = { acidic_flank: 'an acidic flank', basic_flank: 'a basic flank', blocky: 'a blocky IDR', acidic_stretch: 'an acidic stretch',
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
    const views = { overview: overviewView, proteins: proteinsView, flanks: flanksView, idrs: idrsView, sites: sitesView, definitions: definitionsView, about: aboutView };
    const v = view();
    v.innerHTML = '<p class="loading">Loading…</p>';
    C.hideTip();
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

  // ---------------------------------------------------------------- tables (the Methylarginine Site Explorer's)
  // Its markup, styles and helpers (tbl, colHead, colRow, bindRows): click a header to sort (▲ ▼ on the sorted one),
  // Columns to show the optional ones (remembered in this browser), click a row (or Enter on it) to open it, Show more,
  // Download TSV. cols: { k, label, get (value: sorting and TSV), td (cell HTML; default the value), num (right-aligned,
  // sorts high to low first), opt (in the Columns menu; on: shown until the reader changes it), tip, menu }.
  // st: { key, dir, n } is the view's (and its URL's); onState(st) redraws.
  const colsOn = (id, cols) => {
    const s = store.get('cols-' + id);
    if (s) { try { return new Set(JSON.parse(s)); } catch (e) { /* the default */ } }
    return new Set(cols.filter((c) => !c.opt || c.on).map((c) => c.k));
  };
  const PAGE_MORE = 100;
  function tbl(host, o) {
    const { id, cols, rows, st } = o;
    const on = colsOn(id, cols), shown = cols.filter((c) => !c.opt || on.has(c.k));
    const col = cols.find((c) => c.k === st.key) || cols[0];
    const order = rows.slice(), dir = st.dir || (col.num ? -1 : 1);
    const vals = new Map(order.map((i) => [i, col.get(i)]));
    order.sort((a, b) => {
      const x = vals.get(a), y = vals.get(b);
      if (x == null || x === '') return y == null || y === '' ? 0 : 1;
      if (y == null || y === '') return -1;
      return (typeof x === 'string' ? x.localeCompare(y) : x - y) * dir;
    });
    const n = Math.min(order.length, st.n || PAGE_MORE);
    // header cells are upper-cased by the explorer's styles: Greek letters keep their case (δ is not Δ)
    const hc = shown.map((c) => ({ label: esc(c.label).replace(/([α-ω])/g, '<span class="gk">$1</span>'), sort: c.k, cls: c.num ? 'r' : c.cls, tdCls: c.tdCls, tip: c.tip, td: (i) => (c.td ? c.td(i) : esc(c.get(i) ?? '')) }));
    const menu = cols.some((c) => c.opt) ? '<details class="col-menu"><summary class="btn quiet small">Columns</summary><div class="col-pop" role="group" aria-label="Optional columns">' +
      cols.filter((c) => c.opt).map((c) => '<label class="check"><input type="checkbox" data-col="' + c.k + '"' + (on.has(c.k) ? ' checked' : '') + '> ' + esc(c.menu || c.label) + '</label>').join('') +
      '<button class="linkish" type="button" data-col-reset>Default columns</button></div></details>' : '';
    host.innerHTML = '<div class="card tbl-card"><div class="tbl-tools"><span class="muted small">' + (o.note || 'Click a row to open it.') + '</span><span class="tbl-acts">' + menu +
      (o.tsv ? '<button class="btn small" type="button" data-tsv>Download TSV</button>' : '') + '</span></div>' +
      '<div class="tbl-wrap scroll"><table class="tbl stackable" aria-label="' + esc(o.caption) + '">' + C.colHead(hc, { dir }, col.k) + '<tbody>' +
      order.slice(0, n).map((i) => '<tr class="row' + (o.rowCls ? ' ' + o.rowCls(i) : '') + '" tabindex="0" data-i="' + i + '">' + C.colRow(hc, i) + '</tr>').join('') +
      '</tbody></table></div><div class="tbl-foot"><span>Showing ' + fmt(n) + ' of ' + fmt(order.length) + '</span>' +
      (order.length > n ? '<span><button class="btn quiet small" type="button" data-more>Show ' + fmt(Math.min(PAGE_MORE, order.length - n)) + ' more</button>' +
        (order.length <= 3000 ? ' <button class="btn quiet small" type="button" data-all>Show all</button>' : '') + '</span>' : '') + '</div></div>';
    C.bindRows(host, (tr) => o.onOpen(+tr.dataset.i));
    $$('th[data-sort]', host).forEach((th) => th.addEventListener('click', () => {
      const k = th.dataset.sort, c = cols.find((x) => x.k === k);
      o.onState(Object.assign(st, k === col.k ? { dir: -dir } : { key: k, dir: c.num ? -1 : 1 }, { n: PAGE_MORE }));
    }));
    const more = $('[data-more]', host), all = $('[data-all]', host);
    if (more) more.addEventListener('click', () => o.onState(Object.assign(st, { n: n + PAGE_MORE })));
    if (all) all.addEventListener('click', () => o.onState(Object.assign(st, { n: order.length })));
    $$('[data-col]', host).forEach((cb) => cb.addEventListener('change', () => {
      const s = colsOn(id, cols);
      if (cb.checked) s.add(cb.dataset.col); else s.delete(cb.dataset.col);
      store.set('cols-' + id, JSON.stringify([...s]));
      o.onState(st);
      const m = $('.col-menu', host); if (m) m.open = true;
    }));
    const reset = $('[data-col-reset]', host);
    if (reset) reset.addEventListener('click', () => { store.set('cols-' + id, null); o.onState(st); });
    const dl = $('[data-tsv]', host);
    if (dl) dl.addEventListener('click', () => C.download(o.tsv, '# ' + o.caption + (R.scoped ? ' · selection: ' + selectionText() : '') + '\n' + cols.map((c) => c.k).join('\t') + '\n' +
      order.map((i) => cols.map((c) => { const v = c.get(i); return v == null ? '' : String(v); }).join('\t')).join('\n') + '\n'));
  }
  // a table view's sort and length in its URL: sort=key, asc/desc, n
  const tblState = (params, key, dir) => ({ key: params.get('sort') || key, dir: params.has('dir') ? +params.get('dir') : dir, n: +(params.get('n') || PAGE_MORE) });
  const tblParams = (p, st, key, dir) => { if (st.key !== key) p.set('sort', st.key); if (st.dir !== dir) p.set('dir', st.dir); if (st.n > PAGE_MORE) p.set('n', st.n); return p; };

  // ---------------------------------------------------------------- a plot that selects
  // Every row a point on a canvas: all rows faint, the rows that pass the filters coloured by a category, so 100k rows
  // are seen at once. Hover a point: its row (and a ring); click: open it; drag: select a rectangle, which filters the
  // table (and stays in the URL); double-click or Clear: drop it. o: { n, x(i), y(i), rows (passing the filters, box
  // aside), inBox(i), cat(i) -> index into cats, cats: [[label, css colour var]], tip(i) -> html, onOpen(i),
  // box: [x0, x1, y0, y1] | null, onBox(box | null), xLabel, yLabel, xFmt, yFmt, xLog, guides(ctx, X, Y, W, H) }
  // while a plot or map is dragged, no text on the page is selected
  const dragging = (on) => { document.body.classList.toggle('dragging', on); if (on) { try { window.getSelection().removeAllRanges(); } catch (e) { /* none */ } } };
  function scatter(host, o) {
    const W = Math.max(320, Math.round(host.clientWidth || 900)), H = Math.round(Math.min(440, Math.max(300, W * 0.42)));
    const M = { l: 58, r: 16, t: 12, b: 44 }, dpr = window.devicePixelRatio || 1;
    const fx = o.xLog ? (v) => Math.log10(Math.max(1, v)) : (v) => v;
    // domain: all rows, the extremes trimmed (0.1% each side) so one outlier does not squash the rest
    const span = (get, f) => {
      const v = []; for (let i = 0; i < o.n; i++) { const a = get(i); if (a != null && Number.isFinite(a)) v.push(f(a)); }
      v.sort((a, b) => a - b);
      const lo = v[Math.floor(v.length * 0.001)], hi = v[Math.ceil(v.length * 0.999) - 1];
      const pad = (hi - lo) * 0.04 || 1;
      return [o.xMin0 && lo >= 0 ? Math.min(0, lo) : lo - pad, hi + pad];
    };
    const xd = o.xDomain || span(o.x, fx), yd = o.yDomain || span(o.y, (v) => v);
    const X = (v) => M.l + (fx(v) - xd[0]) / (xd[1] - xd[0]) * (W - M.l - M.r), Y = (v) => H - M.b - (v - yd[0]) / (yd[1] - yd[0]) * (H - M.t - M.b);
    const invX = (px) => { const t = xd[0] + (px - M.l) / (W - M.l - M.r) * (xd[1] - xd[0]); return o.xLog ? 10 ** t : t; };
    const invY = (py) => yd[0] + (H - M.b - py) / (H - M.t - M.b) * (yd[1] - yd[0]);
    host.innerHTML = '<div class="sc"><canvas class="sc-base"></canvas><canvas class="sc-over"></canvas></div><div class="sc-key"></div>';
    const [base, over] = $$('canvas', host);
    for (const c of [base, over]) { c.width = W * dpr; c.height = H * dpr; c.style.width = W + 'px'; c.style.height = H + 'px'; }
    const g = base.getContext('2d'), go = over.getContext('2d');
    g.scale(dpr, dpr); go.scale(dpr, dpr);
    const col = (v) => C.cssVar(v) || v, ink = col('--ink-2'), muted = col('--muted'), line = col('--line');
    // axes and grid
    const ticks = (d, log) => {
      if (log) { const out = []; for (let e = Math.ceil(d[0]); e <= Math.floor(d[1]); e++) out.push(e); return out.map((e) => [e, 10 ** e]); }
      const step = [0.01, 0.02, 0.05, 0.1, 0.2, 0.25, 0.5, 1, 2, 5, 10, 20, 25, 50, 100, 200, 500, 1000].find((s) => (d[1] - d[0]) / s <= 7) || 1000;
      const out = []; for (let v = Math.ceil(d[0] / step) * step; v <= d[1] + 1e-9; v += step) out.push([+v.toFixed(6), +v.toFixed(6)]); return out;
    };
    g.font = '11px ' + (C.cssVar('--font-body') || 'sans-serif'); g.lineWidth = 1;
    for (const [t, v] of ticks(xd, o.xLog)) { const px = M.l + (t - xd[0]) / (xd[1] - xd[0]) * (W - M.l - M.r); g.strokeStyle = line; g.beginPath(); g.moveTo(px, M.t); g.lineTo(px, H - M.b); g.stroke(); g.fillStyle = muted; g.textAlign = 'center'; g.fillText((o.xFmt || String)(v), px, H - M.b + 15); }
    for (const [, v] of ticks(yd)) { const py = Y(v); g.strokeStyle = line; g.beginPath(); g.moveTo(M.l, py); g.lineTo(W - M.r, py); g.stroke(); g.fillStyle = muted; g.textAlign = 'right'; g.fillText((o.yFmt || String)(v), M.l - 6, py + 4); }
    g.fillStyle = ink; g.textAlign = 'center'; g.font = '12px ' + (C.cssVar('--font-body') || 'sans-serif');
    g.fillText(o.xLabel, (M.l + W - M.r) / 2, H - 8);
    g.save(); g.translate(14, (M.t + H - M.b) / 2); g.rotate(-Math.PI / 2); g.fillText(o.yLabel, 0, 0); g.restore();
    g.save(); g.beginPath(); g.rect(M.l, M.t, W - M.l - M.r, H - M.t - M.b); g.clip();
    if (o.guides) o.guides(g, X, Y, W, H, M);
    // every row, faint
    g.fillStyle = muted; g.globalAlpha = o.n > 50000 ? 0.12 : 0.2;
    for (let i = 0; i < o.n; i++) { const a = o.x(i), b = o.y(i); if (a == null || b == null) continue; g.fillRect(X(a) - 1, Y(b) - 1, 2, 2); }
    // the rows that pass, in their category's colour; a jitter of a third of a pixel keeps lattice values from one dot
    const shown = o.rows.filter((i) => !o.inBox || o.inBox(i)), big = shown.length > 20000, r = big ? 1.4 : shown.length > 3000 ? 2 : 2.8;
    const cols = o.cats.map((c) => col(c[1])), counts = o.cats.map(() => 0);
    g.globalAlpha = big ? 0.35 : shown.length > 3000 ? 0.55 : 0.8;
    const grid = new Map(), cell = 8;
    for (const i of shown) {
      const a = o.x(i), b = o.y(i);
      if (a == null || b == null) continue;
      const k = o.cat(i), px = X(a), py = Y(b);
      counts[k]++;
      g.fillStyle = cols[k];
      g.beginPath(); g.arc(px, py, r, 0, 6.283); g.fill();
      const key = Math.floor(px / cell) + ',' + Math.floor(py / cell);
      if (!grid.has(key)) grid.set(key, []);
      grid.get(key).push(i);
    }
    g.restore(); g.globalAlpha = 1;
    if (o.box) { const [x0, x1, y0, y1] = o.box; g.strokeStyle = col('--accent'); g.setLineDash([4, 3]); g.strokeRect(X(x0), Y(y1), X(x1) - X(x0), Y(y0) - Y(y1)); g.setLineDash([]); }
    $('.sc-key', host).innerHTML = o.cats.map((c, k) => '<span><i style="background:' + cols[k] + '"></i>' + esc(c[0]) + ' <b>' + fmt(counts[k]) + '</b></span>').join('') +
      '<span class="muted">grey: all ' + fmt(o.n) + ' · drag to select, click a point to open it' + (o.box ? ' · <button class="linkish" type="button" data-clear>clear the selection</button>' : '') + '</span>';
    const clear = $('[data-clear]', host);
    if (clear) clear.addEventListener('click', () => o.onBox(null));
    // hover, click, drag
    const near = (px, py) => {
      let best = null, bd = 100;
      const cx = Math.floor(px / cell), cy = Math.floor(py / cell);
      for (let dx = -1; dx <= 1; dx++) for (let dy = -1; dy <= 1; dy++) for (const i of grid.get((cx + dx) + ',' + (cy + dy)) || []) {
        const d = (X(o.x(i)) - px) ** 2 + (Y(o.y(i)) - py) ** 2;
        if (d < bd) { bd = d; best = i; }
      }
      return best;
    };
    const at = (e) => { const b = over.getBoundingClientRect(); return [e.clientX - b.left, e.clientY - b.top]; };
    let start = null, hover = null;
    const ring = (i) => { go.clearRect(0, 0, W, H); if (i == null) return; go.strokeStyle = col('--ink'); go.lineWidth = 2; go.beginPath(); go.arc(X(o.x(i)), Y(o.y(i)), r + 4, 0, 6.283); go.stroke(); };
    over.addEventListener('pointerdown', (e) => { if (e.button !== 0) return; e.preventDefault(); dragging(true); start = at(e); over.setPointerCapture(e.pointerId); });
    over.addEventListener('pointermove', (e) => {
      const [px, py] = at(e);
      if (start && Math.abs(px - start[0]) + Math.abs(py - start[1]) > 4) {
        C.hideTip(); go.clearRect(0, 0, W, H);
        go.fillStyle = col('--accent'); go.globalAlpha = 0.12; go.fillRect(Math.min(px, start[0]), Math.min(py, start[1]), Math.abs(px - start[0]), Math.abs(py - start[1])); go.globalAlpha = 1;
        go.strokeStyle = col('--accent'); go.strokeRect(Math.min(px, start[0]), Math.min(py, start[1]), Math.abs(px - start[0]), Math.abs(py - start[1]));
        return;
      }
      hover = near(px, py);
      ring(hover);
      over.style.cursor = hover == null ? 'crosshair' : 'pointer';
      if (hover == null) C.hideTip(); else C.showTip(e, o.tip(hover) + '<br><span class="muted">Click to open</span>');
    });
    over.addEventListener('pointerleave', () => { ring(null); C.hideTip(); });
    over.addEventListener('pointerup', (e) => {
      dragging(false);
      if (!start) return;
      const [px, py] = at(e), [sx, sy] = start;
      start = null;
      if (Math.abs(px - sx) + Math.abs(py - sy) > 4) {
        const xs = [invX(Math.min(px, sx)), invX(Math.max(px, sx))], ys = [invY(Math.max(py, sy)), invY(Math.min(py, sy))];
        o.onBox([xs[0], xs[1], ys[0], ys[1]].map((v) => +v.toPrecision(4)));
      } else if (hover != null) { C.hideTip(); o.onOpen(hover); }
    });
    over.addEventListener('dblclick', () => { if (o.box) o.onBox(null); });
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
    const spec = { id: 'ch-kin-' + ctx, title: 'Kinase motifs of phosphosites in ' + ctxLabel(ctx).toLowerCase() + 's', rows, extraCol: true, heads: ['Kinase (motif class)', 'match / sites', '', 'log2 OR'],
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
    const mk = (fl) => W.filter((r) => r.flank === fl).map((r) => ({ key: fl + ':' + r.class, label: label(r.class), sub: fl, value: l2(r.fold), q: r.q, sig: SIG(r.q),
      extra: fmt(r.k) + '/' + fmt(r.M), w: r, title: '<b>' + esc(label(r.class)) + '</b> (' + esc(what(r.class)) + ')<br>' + fmt(r.k) + ' of the ' + fmt(r.n) + ' proteins with an ' + fl + ' flank, against ' + fmt(r.M) + ' of ' + fmt(r.N) +
        ' reviewed proteins<br>fold ' + fx(r.fold) + ', one-sided hypergeometric q ' + fmtP(r.q) }));
    let rows;
    const protSel = R.sel && R.sel.prots;   // an evidence-only selection does not change which proteins are counted
    if (protSel) {
      const res = ST.proteinSets(R.prots);
      rows = res.map((x) => ({ key: x.class, label: label(x.class), value: x.k ? l2(x.fold) : null, q: x.q, sig: SIG(x.q), extra: fmt(x.k) + '/' + fmt(x.M), s: x,
        title: '<b>' + esc(label(x.class)) + '</b><br>' + fmt(x.k) + ' of your ' + fmt(x.n) + ' proteins, against ' + fmt(x.M) + ' of ' + fmt(x.N) + ' reviewed proteins<br>fold ' + (x.fold ? fx(x.fold) : '—') + ', q ' + fmtP(x.q) }));
    } else rows = [{ ph: 'Proteins with an acidic flank' }].concat(mk('acidic'), [{ ph: 'Proteins with a basic flank' }], mk('basic'));
    // section headings ride on the first row of each section, as hbars draws them
    const out = [];
    for (let i = 0; i < rows.length; i++) { if (rows[i].ph) { if (rows[i + 1]) rows[i + 1].head = rows[i].ph; continue; } out.push(rows[i]); }
    for (const r of out) r.main = { or: r.value == null ? null : 2 ** r.value, q: r.q, sig: r.sig, fold: true };
    return { id: 'ch-classes', title: 'Protein classes with an acidic or basic flank', rows: out, fold: true, extraCol: true, heads: ['Class', 'in set / all', '', 'log2 fold'],
      subHtml: (protSel ? 'Your proteins' : 'Proteins with a flank in the top tenth of the acidic or the basic score') + ' against all ' + fmt(B.counts.proteins) + ' reviewed human proteins: one-sided hypergeometric test, fold enrichment (k/n ÷ M/N), BH across the classes.' + docLink('tests'),
      subText: (protSel ? 'Your proteins' : 'Proteins with an acidic or a basic flank') + ' against all reviewed human proteins: one-sided hypergeometric, fold enrichment, BH.',
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
    const spec = { id: 'ch-charge', title: 'Charge blocks after modification, against random placement', rows, extraCol: true, heads: ['Modification', 'sharper / IDRs', '', 'log2 odds'],
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
    return { id: 'ch-crosstalk', title: 'Methylarginines by offset from a phosphosite', rows, keepOrder: true, extraCol: true, heads: ['Arginine at', 'methyl-R', '', 'log2 OR'],
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
    return { id: 'ch-partners', title: 'Partner pairs by flank class', rows, keepOrder: true, extraCol: true, heads: ['Partner pair', 'edges', '', 'log2 O/E'],
      subHtml: 'Partner pairs by the flank class of their proteins (acidic, basic, both, neither), observed against 2,000 shuffles of the classes within degree decile and nuclear status; histone-like proteins left out. The contrast compares opposite-charge pairs with same-charge pairs. Two-sided permutation p. Whole dataset only.' + docLink('tests'),
      subText: 'Partner pairs by flank class, observed over expected under class shuffles within degree and nuclear status; histones left out.',
      axis: 'log2 observed / expected partner pairs', valueLabel: 'log2_obs_over_exp' };
  }
  function stoichSpec(which) {
    const label = which === 'occ' ? 'Phosphosite occupancy (Koenig et al. 2026)' : 'Acetyl-K stoichiometry (Hansen et al. 2019)';
    const measure = which === 'occ' ? 'phospho-S/T occupancy (Koenig 2026)' : 'acetyl-K stoichiometry (Hansen 2019)';
    const ctxs = ['acidic_flank', 'basic_flank', 'flank', 'idr', 'blocky', 'basic_patch', 'domain'];
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
    return { id: 'ch-agnostic-' + measure, title: label + ' by depth from the IDR end, over the same IDR shuffled', rows: [], figSvg: () => agSvg(measure, false), figCss: AG_CSS,
      subHtml: 'Every IDR read from each end, each depth against the same IDR’s residues in random order (20 shuffles), so composition and length are kept and only position is removed. Lines: the IDR ends that touch each kind of domain, and the ends at the protein’s own terminus. Filled points: 95% CI (bootstrap over proteins) excludes 0. Whole dataset only.' + docLink('agnostic'),
      subText: label + ' (' + unit + ') by distance from the IDR end, against the same IDR shuffled; 95% CI by bootstrap over proteins.',
      figTsv: () => '# ' + label + ', ' + unit + '; excess over the same IDR shuffled\nend\tbin\texcess\tci_low\tci_high\tends\tproteins\n' + rows.map((r) => [r.end, r.bin, r.excess * scale, r.ci_low * scale, r.ci_high * scale, r.ends, r.proteins].join('\t')).join('\n') + '\n' };
  }

  // ---------------------------------------------------------------- which domain families carry a charged extension? (dtidr_family_profiles.py)
  // Every IDR end touching an InterPro entry, read from the end, each depth band against a reference (the same IDR, or
  // the protein's non-domain residues); per sequence cluster (paralogs collapsed), t-test and BH, length from the bands.
  const FAM_REF = [['idr', 'against its own IDR'], ['protein', 'against the protein’s non-domain residues']];
  const FAM_MEASURES = [['basic', 'K/R', 'basic'], ['acidic', 'D/E', 'acidic'], ['net_charge', 'net charge', 'net']];
  // the tile and notes read overview.json (either side); the section reads families.json (every side, band and measure)
  const famRows = (ref, side) => tableRows(FAM && FAM.families || OV.families).filter((r) => r.reference === ref && r.side === side);
  const famSig = (r) => FAM_MEASURES.some(([m]) => SIG(r[m + '_proximal_q']));
  const pts = (v, d = 1) => (v == null ? '—' : (v >= 0 ? '+' : '−') + Math.abs(100 * v).toFixed(d));
  function familySpec(ref, side, measure) {
    const [, mlab] = FAM_MEASURES.find((m) => m[0] === measure);
    const all = famRows(ref, side);
    const rows = all.filter((r) => r[measure + '_proximal'] != null).sort((a, b) => b[measure + '_proximal'] - a[measure + '_proximal']).map((r) => ({
      key: r.interpro, label: r.name, value: 100 * r[measure + '_proximal'], q: r[measure + '_proximal_q'], sig: SIG(r[measure + '_proximal_q']), extra: fmt(r.units), w: r,
      title: '<b>' + esc(r.name) + '</b> (' + esc(r.interpro) + ', ' + esc(r.type) + ')<br>' + mlab + ' at residues 1–10 of the IDR ends beside it: ' + pts(r[measure + '_proximal']) + ' points ' + esc(FAM_REF.find((x) => x[0] === ref)[1]) +
        ', q ' + fmtP(r[measure + '_proximal_q']) + '<br>11–30: ' + pts(r[measure + '_mid']) + ' points, q ' + fmtP(r[measure + '_mid_q']) + (r[measure + '_length'] ? '<br>extension length about ' + r[measure + '_length'] + ' residues (bands whose 95% CI excludes 0)' : '') +
        '<br>' + fmt(r.ends) + ' ends on ' + fmt(r.proteins) + ' proteins, ' + fmt(r.units) + ' after collapsing paralogs' }));
    for (const r of rows) r.main = { or: 2 ** (r.value / 10), q: r.q, sig: r.sig };
    return { id: 'ch-fam-' + ref + '-' + side + '-' + measure, title: mlab + ' excess at residues 1–10 beside each domain family', rows, extraCol: true, heads: ['InterPro entry', 'units', '', 'points'],
      format: (v) => (v == null ? '—' : (v >= 0 ? '+' : '−') + Math.abs(v).toFixed(1)), tick: (t) => (t > 0 ? '+' : '') + t,
      subHtml: 'Each InterPro domain or repeat with IDR ends on at least 10 independent sequences: the ' + esc(mlab) + ' fraction of residues 1–10 of the IDR, read from its end at the domain, minus the same fraction ' + esc(FAM_REF.find((x) => x[0] === ref)[1]) +
        ' (percentage points). One-sample t-test over sequence clusters (paralogs sharing a 12-residue stretch count once), BH across entries, sides, measures and bands. Whole dataset only.' + docLink('families'),
      subText: mlab + ' fraction at residues 1-10 of IDR ends beside each InterPro entry, minus the reference (' + FAM_REF.find((x) => x[0] === ref)[1] + '); t-test over sequence clusters, BH.',
      axis: 'percentage points of ' + mlab + ' at residues 1–10, ' + FAM_REF.find((x) => x[0] === ref)[1], valueLabel: 'excess_points',
      tsvCols: [['interpro', (r) => r.key], ['q', (r) => r.q], ['mid_points', (r) => 100 * r.w[measure + '_mid']], ['mid_q', (r) => r.w[measure + '_mid_q']], ['length', (r) => r.w[measure + '_length']], ['proteins', (r) => r.w.proteins], ['units', (r) => r.w.units]] };
  }
  // one family's bands: K/R and D/E excess by depth with 95% CIs (bootstrap over sequence clusters)
  function famBandSvg(ipr, ref, side) {
    const rows = tableRows(FAM.family_bands).filter((r) => r.interpro === ipr && r.reference === ref && r.side === side);
    const bands = [...new Set(rows.map((r) => r.band))], W = 560, H = 250, L = 56, Rm = 120, T = 12, Bm = 42;
    const vals = rows.filter((r) => r.measure !== 'net_charge').flatMap((r) => [r.ci_low, r.ci_high]).filter((v) => v != null).map((v) => 100 * v).concat([0]);
    let lo = Math.min(...vals), hi = Math.max(...vals); const pad = (hi - lo) * 0.08 || 1; lo -= pad; hi += pad;
    const X = (i) => L + (i + 0.5) / bands.length * (W - L - Rm), Y = (v) => T + (1 - (v - lo) / (hi - lo)) * (H - T - Bm);
    const step = [1, 2, 5, 10, 20].find((s) => (hi - lo) / s <= 6) || 25;
    let g = '';
    for (let v = Math.ceil(lo / step) * step; v <= hi; v += step) g += '<line class="ag-grid" x1="' + L + '" x2="' + (W - Rm) + '" y1="' + Y(v).toFixed(1) + '" y2="' + Y(v).toFixed(1) + '"/><text class="ag-tick" x="' + (L - 6) + '" y="' + (Y(v) + 4).toFixed(1) + '" text-anchor="end">' + (v > 0 ? '+' : '') + v + '</text>';
    bands.forEach((b, i) => { g += '<text class="ag-tick" x="' + X(i).toFixed(1) + '" y="' + (H - Bm + 16) + '" text-anchor="middle">' + b + '</text>'; });
    g += '<text class="ag-axt" x="' + ((L + W - Rm) / 2) + '" y="' + (H - 6) + '" text-anchor="middle">residues from the IDR’s end at the domain</text><text class="ag-axt" transform="translate(14 ' + ((T + H - Bm) / 2) + ') rotate(-90)" text-anchor="middle">excess, percentage points</text>';
    g += '<line class="ag-zero" x1="' + L + '" x2="' + (W - Rm) + '" y1="' + Y(0).toFixed(1) + '" y2="' + Y(0).toFixed(1) + '"/>';
    const name = (tableRows(FAM.families).find((r) => r.interpro === ipr) || {}).name || ipr;
    for (const [m, lab, cls, dx] of [['basic', 'K/R', 'fb', -5], ['acidic', 'D/E', 'fa', 5]]) {
      const pts2 = bands.map((b, i) => { const r = rows.find((x) => x.band === b && x.measure === m); return r && r.excess != null ? [X(i) + dx, r] : null; }).filter(Boolean);
      if (!pts2.length) continue;
      g += '<polyline class="ag-line fam-' + cls + '" points="' + pts2.map(([x, r]) => x.toFixed(1) + ',' + Y(100 * r.excess).toFixed(1)).join(' ') + '"/>';
      for (const [x, r] of pts2) {
        const sig = r.ci_low > 0 || r.ci_high < 0;
        g += '<line class="ag-ci fam-' + cls + '" x1="' + x.toFixed(1) + '" x2="' + x.toFixed(1) + '" y1="' + Y(100 * r.ci_low).toFixed(1) + '" y2="' + Y(100 * r.ci_high).toFixed(1) + '"/>' +
          '<circle class="ag-pt fam-' + cls + (sig ? '' : ' ns') + '" cx="' + x.toFixed(1) + '" cy="' + Y(100 * r.excess).toFixed(1) + '" r="4" data-tip="' + esc('<b>' + esc(name) + ', residues ' + r.band + '</b><br>' + lab + ' ' + pts(r.excess) + ' points, 95% CI ' + pts(r.ci_low) + ' to ' + pts(r.ci_high) + ' · ' + fmt(r.units) + ' sequence clusters') + '"/>';
      }
      const last = pts2[pts2.length - 1];
      g += '<text class="ag-lab fam-' + cls + '-t" x="' + (W - Rm + 8) + '" y="' + (Y(100 * last[1].excess) + 4).toFixed(1) + '">' + lab + '</text>';
    }
    return '<svg class="ag-plot" viewBox="0 0 ' + W + ' ' + H + '" role="img" aria-label="' + esc('K/R and D/E excess by depth beside ' + name) + '">' + g + '</svg>';
  }
  const FAM_CSS = (pal) => AG_CSS(pal) + '.fam-fb{stroke:#2F6FD6;fill:#2F6FD6}.fam-fa{stroke:#D64545;fill:#D64545}.ag-line.fam-fb,.ag-line.fam-fa{fill:none}.fam-fb-t{fill:#2F6FD6}.fam-fa-t{fill:#D64545}';
  function famMini() {
    const rows = famRows('idr', 'both').filter((r) => SIG(r.basic_proximal_q) || SIG(r.acidic_proximal_q)).sort((a, b) => b.basic_proximal - a.basic_proximal).slice(0, 6);
    const mx = Math.max(...rows.map((r) => Math.abs(r.basic_proximal))) || 1;
    return '<div class="fam-mini">' + rows.map((r) => '<div class="fm-row" data-tip="' + esc('<b>' + esc(r.name) + '</b><br>K/R ' + pts(r.basic_proximal) + ' points at residues 1–10 over the rest of the IDR, q ' + fmtP(r.basic_proximal_q) + (r.basic_length ? '; about ' + r.basic_length + ' residues long' : '') + '<br>D/E ' + pts(r.acidic_proximal) + ' points, q ' + fmtP(r.acidic_proximal_q)) + '">' +
      '<span class="fm-l">' + esc(r.name.replace(/ domain$/, '')) + '</span><span class="fm-t"><i class="' + (r.basic_proximal > 0 ? 'b' : 'a') + '" style="width:' + (100 * Math.abs(r.basic_proximal) / mx).toFixed(1) + '%"></i></span><span class="fm-v">' + pts(r.basic_proximal) + '</span></div>').join('') + '</div>';
  }

  // ---------------------------------------------------------------- how are IDRs and flanks charged? (idrkit.charge_type)
  // Groups of IDRs by what they touch, and flank windows; recomputed from idrs.json for a selection of proteins.
  const CT_GROUP = { all: 'All IDRs', acidic: 'IDRs at an acidic flank', basic: 'IDRs at a basic flank', flank: 'IDRs at any domain flank', none: 'IDRs away from domain flanks',
    fl_all: 'All flank windows', fl_acidic: 'Acidic flank windows', fl_basic: 'Basic flank windows' };
  const CT_SHOW = [3, 4, 2, 1];   // acidic, basic, mixed, few: the order of every stack
  function ctGroups() {
    const prots = R.sel && R.sel.prots;
    if (!prots || !ID) return { groups: OV.charge_types.idr_groups.map(([k, n]) => ({ key: k, label: CT_GROUP[k], n })), whole: true };
    const cnt = Object.fromEntries(Object.keys(CT_GROUP).slice(0, 5).map((k) => [k, [0, 0, 0, 0]]));
    for (let i = 0; i < ID.p.length; i++) {
      if (!prots.has(ID.p[i])) continue;
      for (const k of ['all', ID.top[i] & 1 ? 'acidic' : null, ID.top[i] & 2 ? 'basic' : null, ID.zb[i] != null ? 'flank' : 'none'].filter(Boolean)) cnt[k][ID.ct[i] - 1]++;
    }
    return { groups: Object.keys(cnt).map((k) => ({ key: k, label: CT_GROUP[k], n: cnt[k] })), whole: false };
  }
  // one stacked bar per group, with direct percentages: by default the charge words (acidic, basic, mixed, few), or any
  // categories (o.cats: [[index into n, label, css class, definition]])
  function ctStackSvg(groups, o = {}) {
    const cats = o.cats || CT_SHOW.map((t) => [t - 1, CT[t], 'cs' + t, CT_DEF[t]]);
    const W = o.W || 780, ML = o.ML || 240, RH = o.compact ? 16 : 26, T = 4, bw = W - ML - 8;
    // the key wraps onto as many lines as the width needs
    const keyAt = [];
    let lx = ML, ly = 0;
    for (const [, lab] of o.compact ? [] : cats) { const w = 30 + C.textW(lab, 11); if (lx > ML && lx + w > W) { lx = ML; ly++; } keyAt.push([lx, ly]); lx += w; }
    const H = T + groups.length * RH + (o.compact ? 2 : 12 + (ly + 1) * 16);
    let g = '';
    groups.forEach((gr, k) => {
      const tot = gr.n.reduce((a, b) => a + b, 0) || 1, y = T + k * RH;
      g += '<text class="ct-lab" x="' + (ML - 8) + '" y="' + (y + RH / 2 + 4) + '" text-anchor="end">' + esc(gr.label) + (o.compact ? '' : ' (' + fmt(tot) + ')') + '</text>';
      let x = ML;
      for (const [j, lab, cls, def] of cats) {
        const c = gr.n[j];
        if (!c) continue;
        const w = bw * c / tot;
        g += '<rect class="ct-seg ' + cls + '" x="' + x.toFixed(1) + '" y="' + (y + 3) + '" width="' + Math.max(0.5, w - 1.5).toFixed(1) + '" height="' + (RH - 6) + '" rx="2" data-tip="' + esc('<b>' + esc(gr.label) + '</b><br>' + esc(lab) + ': ' + fmt(c) + ' of ' + fmt(tot) + ' (' + (100 * c / tot).toFixed(1) + '%)<br>' + esc(def)) + '"/>';
        if (!o.compact && w > 34) g += '<text class="ct-in inv" x="' + (x + 4).toFixed(1) + '" y="' + (y + RH / 2 + 4) + '">' + Math.round(100 * c / tot) + '%</text>';
        x += w;
      }
    });
    cats.forEach(([, lab, cls], k) => { if (o.compact) return; const [x, l] = keyAt[k], y = T + groups.length * RH + 10 + l * 16; g += '<rect class="ct-seg ' + cls + '" x="' + x.toFixed(1) + '" y="' + y + '" width="12" height="10" rx="2"/><text class="ct-lab" x="' + (x + 16).toFixed(1) + '" y="' + (y + 9) + '">' + esc(lab) + '</text>'; });
    return '<svg class="ct-plot' + (o.compact ? ' compact' : '') + '" viewBox="0 0 ' + W + ' ' + H + '" role="img" aria-label="' + esc(o.aria || 'Charge of IDRs and flanks by group') + '">' + g + '</svg>';
  }
  const CT_CSS = (pal) => '.ct-lab{font-size:11px;fill:' + pal.ink + '}.ct-in{font-size:10px;fill:' + pal.ink + '}.ct-in.inv{fill:#fff}.cs1{fill:#8A979F}.cs2{fill:' + TRACK_LIGHT.hybrid + '}.cs3{fill:' + TRACK_LIGHT.acidic + '}.cs4{fill:' + TRACK_LIGHT.basic + '}' +
    [['#8A979F'], ['#C98500'], ['#8E44AD'], ['#D64545'], ['#2F6FD6']].map(([c], k) => '.ds' + (k + 1) + '{fill:' + c + '}').join('');
  // the same groups in Das & Pappu's classes (localCIDER), the field's standard, with their predicted shapes
  function dpSpec() {
    const DPG = { all: 'All IDRs', acidic: 'IDRs at an acidic flank', basic: 'IDRs at a basic flank', flank: 'IDRs at any domain flank', none: 'IDRs away from domain flanks', fl_all: 'All flank windows', fl_acidic: 'Acidic flank windows', fl_basic: 'Basic flank windows' };
    const prots = R.sel && R.sel.prots;
    let groups = OV.cider.idr_groups.map((g) => ({ key: g[0], label: DPG[g[0]], n: g[2] }));
    if (prots && ID) {
      const cnt = Object.fromEntries(['all', 'acidic', 'basic', 'flank', 'none'].map((k) => [k, [0, 0, 0, 0, 0]]));
      for (let i = 0; i < ID.p.length; i++) if (prots.has(ID.p[i])) for (const k of ['all', ID.top[i] & 1 ? 'acidic' : null, ID.top[i] & 2 ? 'basic' : null, ID.zb[i] != null ? 'flank' : 'none'].filter(Boolean)) cnt[k][(ID.reg[i] || 1) - 1]++;
      groups = Object.keys(cnt).map((k) => ({ key: k, label: DPG[k], n: cnt[k] }));
    } else groups = groups.concat(OV.cider.flank_groups.map((g) => ({ key: g[0], label: DPG[g[0]], n: g[2] })));
    const cats = [1, 2, 3, 4, 5].map((r) => [r - 1, REGION_SHORT[r], 'ds' + r, REGION[r] + ': ' + REGION_RULE[r] + '; predicted shape ' + REGION_SHAPE[r]]);
    return { id: 'ch-dp', title: 'Das–Pappu classes of IDRs and flank windows', rows: [], figSvg: () => ctStackSvg(groups, { cats, aria: 'Das–Pappu classes by group' }), figCss: CT_CSS,
      subHtml: 'The same groups in the five classes of Das & Pappu’s diagram of states (localCIDER), named in plain words: R1 weakly charged chains are compact, R4 and R5 strongly charged ones of one sign expand, R3 mixed ones fold back or not depending on how their charges are mixed (κ), R2 lies between. Hover a segment for the rule.' + docLink('cider'),
      subText: 'Das-Pappu diagram-of-states class (localCIDER) of each IDR and flank window: R1 FCR<0.25, R2 0.25-0.35, R3 FCR>0.35 and |NCPR|<0.35, R4 f->0.35, R5 f+>0.35.',
      figTsv: () => '# Das-Pappu classes' + (prots ? ' · selection: ' + selectionText() : '') + '\ngroup\tR1\tR2\tR3\tR4\tR5\n' + groups.map((g) => [g.label].concat(g.n).join('\t')).join('\n') + '\n' };
  }
  function ctSpec() {
    const cg = ctGroups(), fl = OV.charge_types.flank_groups.map(([k, n]) => ({ key: k, label: CT_GROUP[k], n }));
    const show = cg.groups.concat(cg.whole ? fl : []);
    return { id: 'ch-ctype', title: 'Charge of IDRs and flank windows', rows: [], figSvg: () => ctStackSvg(show), figCss: CT_CSS, meta: cg,
      subHtml: 'Each IDR (and each 50-residue flank window) in words: <b>acidic</b> at a net charge of −0.10 per residue or below, <b>basic</b> at +0.10 or above; otherwise <b>mixed</b> when a quarter or more of its residues are charged, <b>few charges</b> when fewer.' + (cg.whole ? '' : ' Your proteins’ IDRs; flank windows for the whole dataset only.') + docLink('ctype'),
      subText: 'Charge type of each IDR and flank window: acidic (net charge per residue <= -0.10), basic (>= +0.10), mixed (fraction charged >= 0.25), few charges.',
      figTsv: () => '# Charge types' + (cg.whole ? '' : ' · selection: ' + selectionText()) + '\ngroup\tacidic\tbasic\tmixed\tfew\n' + show.map((g) => [g.label].concat(CT_SHOW.map((t) => g.n[t - 1])).join('\t')).join('\n') + '\n' };
  }
  // per type: IDRs at an acidic (basic) flank against the other IDRs, pooled Fisher, BH across the types
  function ctTestSpec() {
    const cg = ctGroups(), get = (k) => cg.groups.find((g) => g.key === k).n, all = get('all');
    const rows = [];
    for (const [k, head] of [['acidic', 'IDRs at an acidic flank'], ['basic', 'IDRs at a basic flank']]) {
      const inn = get(k), ti = inn.reduce((a, b) => a + b, 0), ta = all.reduce((a, b) => a + b, 0), sub = [];
      for (const t of CT_SHOW) {
        const a = inn[t - 1], b = ti - a, c = all[t - 1] - a, d = ta - ti - c, oc = C.oddsCi(a, b, c, d);
        sub.push({ key: k + t, head: sub.length ? null : head, label: CT[t], value: oc.or ? +Math.log2(oc.or).toFixed(3) : null, p: C.fisherP(a, b, c, d), extra: fmt(a) + '/' + fmt(ti), oc,
          title: '<b>' + esc(head) + ': ' + CT[t].toLowerCase() + '</b><br>' + fmt(a) + ' of ' + fmt(ti) + ' (' + (100 * a / (ti || 1)).toFixed(1) + '%) against ' + fmt(c) + ' of ' + fmt(ta - ti) + ' other IDRs (' + (100 * c / ((ta - ti) || 1)).toFixed(1) + '%)<br>OR ' + C.fmtOr(oc.or, oc.lo, oc.hi) });
      }
      C.bhAdjust(sub, (x) => x.p);
      for (const x of sub) { x.sig = SIG(x.q); x.title += ', Fisher q ' + fmtP(x.q); x.main = { or: x.oc.or, lo: x.oc.lo, hi: x.oc.hi, q: x.q, sig: x.sig }; }
      rows.push(...sub);
    }
    return { id: 'ch-ctype-test', title: 'Charge of the IDRs at acidic and basic flanks', rows, keepOrder: true, extraCol: true, heads: ['Charge', 'IDRs', '', 'log2 OR'],
      subHtml: 'IDRs overlapping an acidic (basic) flank against all other IDRs' + (cg.whole ? '' : ' of your proteins') + ': the odds of each charge type. Pooled Fisher test (IDRs are the unit), BH across the types. The flank scores read charge from the same residues, so acidic IDRs at acidic flanks are expected; the size of the rest is the information.' + docLink('ctype'),
      subText: 'IDRs overlapping an acidic (basic) flank against all other IDRs: odds of each charge type, pooled Fisher, BH.', axis: 'log2 odds ratio vs other IDRs', valueLabel: 'log2_odds_ratio' };
  }

  // ---------------------------------------------------------------- Overview
  // the tiles: key, group, question (the section's heading); titles are nouns (TILE_TITLE)
  const TILES = [
    ['agnostic', 'Is a domain flank a unit?', 'Is the sequence beside a folded domain different from the rest of its IDR?'],
    ['families', 'Is a domain flank a unit?', 'Which domain families have a charged extension, and how long is it?'],
    ['ctype', 'Is a domain flank a unit?', 'How charged are IDRs and the flanks of domains?'],
    ['acidic_flank', 'Domain flanks', 'Which modifications sit in acidic flanks?'], ['basic_flank', 'Domain flanks', 'Which modifications sit in basic flanks?'], ['flank', 'Domain flanks', 'Which modifications sit beside a folded domain?'],
    ['blocky', 'Charge patterning', 'Which modifications sit in IDRs with charge blocks?'], ['charge', 'Charge patterning', 'Do modifications sharpen or soften charge blocks?'],
    ['acidic_stretch', 'Sequence elements', 'Which modifications sit in acidic stretches?'], ['basic_patch', 'Sequence elements', 'Which modifications sit in basic patches?'], ['rg_cluster', 'Sequence elements', 'Which modifications sit in RG/RGG clusters?'],
    ['basophilic', 'Sequence elements', 'Are the arginines and lysines just before a phosphosite modified?'],
    ['idr', 'Order and disorder', 'Which modifications sit in IDRs?'], ['domain', 'Order and disorder', 'Which modifications sit in folded domains?'],
    ['kin', 'Writers and occupancy', 'Which kinases’ motifs do the phosphosites of each flank match?'], ['occupancy', 'Writers and occupancy', 'Are sites in acidic or basic flanks more fully modified?'],
    ['classes', 'Proteins and partners', 'What kinds of protein have acidic or basic flanks?'], ['crosstalk', 'Proteins and partners', 'Are arginines next to phosphosites methylated?'], ['partners', 'Proteins and partners', 'Do acidic- and basic-flank proteins bind each other?'],
  ];
  const TILE_TITLE = { agnostic: 'Flank profiles', families: 'Domain families', ctype: 'Charge of IDRs and flanks', kin: 'Kinase motifs', charge: 'Charge blocks after modification', occupancy: 'Occupancy', classes: 'Protein classes', crosstalk: 'Phospho–methyl crosstalk', partners: 'Partners' };
  const tileTitle = (key) => (CTX_ORDER().includes(key) ? ctxLabel(key) + (key.endsWith('flank') || key === 'flank' ? 's' : '') : TILE_TITLE[key]);
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
      const a = kinaseSpec('acidic_flank'), b = kinaseSpec('basic_flank');
      return { spec: kinaseSpec(UI.kinCtx || 'acidic_flank'), groups: [['Acidic flanks', a.rows], ['Basic flanks', b.rows]], rows: a.rows.concat(b.rows), specs: [a, b] };
    }
    if (key === 'charge') { const spec = chargeSpec(); return { spec, rows: spec.rows }; }
    if (key === 'classes') { const spec = classSpec(); return { spec, rows: spec.rows, fold: true }; }
    if (key === 'crosstalk') { const spec = crosstalkSpec(); return { spec, rows: spec.rows, whole: true }; }
    if (key === 'partners') { const spec = partnersSpec(); return { spec, rows: spec.rows, whole: true }; }
    if (key === 'occupancy') { const occ = stoichSpec('occ'), ac = stoichSpec('ac'); return { spec: occ, spec2: ac, occupancy: true }; }
    if (key === 'agnostic') return { spec: agSpec(UI.agMeasure || 'basic'), agnostic: true, whole: true };
    if (key === 'families') return { spec: familySpec(UI.famRef || 'idr', UI.famSide || 'both', UI.famMeasure || 'basic'), families: true, whole: true };
    if (key === 'ctype') { const spec = ctSpec(); return { spec, spec2: dpSpec(), spec3: ctTestSpec(), ctype: true, rows: [] }; }
    return null;
  }
  function tileOf(key, sec) {
    if (sec.agnostic) {
      const r = tableRows(OV.agnostic), get = (end, bin) => r.find((x) => x.measure === 'basic' && x.end === end && x.bin === bin) || {};
      const dna = get('domain:dna', '1-5'), dna2 = get('domain:dna', '16-20'), oth = get('domain:other', '1-5');
      return { vis: agSvg('basic', true) + '<div class="ag-key">' + AG_SERIES.map(([, n, c]) => '<span><i class="ag-sw ag-' + c + '"></i>' + esc(n) + '</span>').join('') + '</div>',
        head: 'Only beside DNA-binding domains: <b>+' + (100 * dna.excess).toFixed(1) + '</b> points of K/R at residues 1–5, gone by 16–20 (' + (100 * dna2.excess).toFixed(1) + '); other domains +' + (100 * oth.excess).toFixed(1) + ' at 1–5 only' };
    }
    if (sec.families) {
      const both = famRows('idr', 'both'), up = both.filter((r) => SIG(r.basic_proximal_q) && r.basic_proximal > 0), ac = both.filter((r) => SIG(r.acidic_proximal_q) && r.acidic_proximal > 0);
      return { vis: famMini(), head: fmt(up.length) + ' of ' + fmt(both.length) + ' families carry a basic extension (DNA-binding folds first); ' + (ac.length ? fmt(ac.length) + ' an acidic one' : '<b>none</b> an acidic one') };
    }
    if (sec.ctype) {
      const G = sec.spec.meta.groups, sh = (k, t) => { const g = G.find((x) => x.key === k); const n = g.n.reduce((a, b) => a + b, 0) || 1; return Math.round(100 * g.n[t - 1] / n); };
      const short = { all: 'All IDRs', acidic: 'At acidic flanks', basic: 'At basic flanks' };
      return { vis: ctStackSvg(G.filter((g) => short[g.key]).map((g) => Object.assign({}, g, { label: short[g.key] })), { compact: true, W: 300, ML: 92 }),
        head: 'Acidic: <b>' + sh('acidic', 3) + '%</b> of IDRs at acidic flanks, ' + sh('all', 3) + '% of all IDRs; basic: ' + sh('basic', 4) + '% at basic flanks, ' + sh('all', 4) + '% of all' };
    }
    if (sec.occupancy) {
      const rows = sec.spec.meta, mx = Math.max(...rows.map((r) => (r.x.median_in || 0)), ...rows.map((r) => r.x.median_out || 0));
      const meters = rows.filter((r) => ['acidic_flank', 'basic_flank', 'idr', 'domain'].includes(r.c)).map((r) => ({ label: ctxLabel(r.c), v: (r.x.median_in || 0) / mx, ref: (r.x.median_out || 0) / mx,
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
      (R.scoped ? ' · <b>' + fmt(R.sites.length) + '</b> sites on <b>' + fmt(R.prots.length) + '</b> proteins' : ' · every test and table runs again on a selection') + '</summary>' +
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
    if (sub === 'families') { await loadFam(); if (token !== routeToken) return; }
    if (sub === 'ctype' && R.sel && R.sel.prots) await loadIdrs();
    if (sub) { sectionView(sub); return; }
    C.A_ROWS.length = 0;
    const ov = {};
    for (const [key, group, q] of TILES) {
      const sec = buildSection(key);
      if (!sec) continue;
      SECTIONS[key] = sec;
      C.foldHtml(key, tileTitle(key), '', false, '', null, { q, how: sec.whole ? 'whole' : 'retest' });
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
      [c.flanks, 'domain flanks (' + fmt(c.acidic_flanks) + ' acidic, ' + fmt(c.basic_flanks) + ' basic)', 'flanks'], [c.idrs, 'IDRs (≥ ' + OV.settings.min_idr + ' aa)', 'idrs'], [c.proteins, 'reviewed proteins', 'proteins'], [c.sites, 'modification sites', 'proteins'],
      [c.by_kind.ph, 'phospho-S/T'], [c.by_kind.mR + c.by_kind.cit, 'methylated or citrullinated R'], [c.by_kind.kac + c.by_kind.ub + c.by_kind.sumo, 'acetyl, ubiquitin, SUMO-K'], [c.by_kind.adpr + c.by_kind.glu, 'ADP-ribosyl and glutamyl sites'],
    ].map(([n, l, go]) => (go ? '<a class="stat" href="#' + go + '">' : '<div class="stat">') + '<b>' + fmt(n) + '</b><span>' + esc(l) + '</span>' + (go ? '</a>' : '</div>')).join('') + '</div>' + C.aHome(lead, ov) + ideasHtml();
    $$('.a-tile[data-go], .a-hi-c[data-go]').forEach((b) => b.addEventListener('click', () => { const s = selParams().toString(); location.hash = 'overview/' + b.dataset.go + (s ? '?' + s : ''); }));
    bindSelPanel('overview');
  }

  function sectionView(key) {
    const sec = buildSection(key);
    if (!sec) { view().innerHTML = '<div class="panel">No section called ' + esc(key) + '. <a href="#overview">Back to the Overview</a></div>'; return; }
    const idx = TILES.findIndex((t) => t[0] === key), t = TILES[idx];
    const s = selParams().toString(), q = s ? '?' + s : '';
    const nav = '<nav class="sec-nav"><a href="#overview' + q + '">‹ All tiles</a>' + (idx > 0 ? '<a href="#overview/' + TILES[idx - 1][0] + q + '">‹ ' + esc(tileTitle(TILES[idx - 1][0])) + '</a>' : '') +
      (idx < TILES.length - 1 ? '<a href="#overview/' + TILES[idx + 1][0] + q + '">' + esc(tileTitle(TILES[idx + 1][0])) + ' ›</a>' : '') + '</nav>';
    const sel = (id, label, opts, cur) => '<label>' + esc(label) + ' <select id="' + id + '">' + opts.map(([k, l]) => '<option value="' + k + '"' + (k === cur ? ' selected' : '') + '>' + esc(l) + '</option>').join('') + '</select></label>';
    let body = '';
    if (key === 'kin') {
      sec.spec.tools = sel('kinctx', 'Context', KIN_CTX.map((c) => [c, ctxLabel(c)]), UI.kinCtx || 'acidic_flank');
      body = card(sec.spec);
    } else if (sec.occupancy) body = card(sec.spec) + card(sec.spec2);
    else if (sec.ctype) body = card(sec.spec) + card(sec.spec2) + card(sec.spec3);
    else if (sec.agnostic) {
      sec.spec.tools = sel('agm', 'Measure', AG_MEASURES.map(([k, l]) => [k, l]), UI.agMeasure || 'basic');
      body = card(sec.spec, '<div class="ag-key">' + AG_SERIES.map(([, n, c]) => '<span><i class="ag-sw ag-' + c + '"></i>' + esc(n) + '</span>').join('') + '</div>') + agSummaryHtml();
    } else if (sec.families) {
      const ref = UI.famRef || 'idr', side = UI.famSide || 'both', m = UI.famMeasure || 'basic';
      sec.spec.tools = sel('famref', 'Reference', FAM_REF, ref) + sel('famside', 'Side', [['both', 'either side'], ['N', 'IDR before the domain'], ['C', 'IDR after the domain']], side) + sel('famm', 'Measure', FAM_MEASURES.map(([k, l]) => [k, l]), m);
      const pick = UI.famPick && famRows(ref, side).some((r) => r.interpro === UI.famPick) ? UI.famPick : (sec.spec.rows[0] || {}).key;
      const pr = famRows(ref, side).find((r) => r.interpro === pick) || {};
      const band = { id: 'ch-fam-bands', title: (pr.name || '') + ': K/R and D/E by depth', rows: [], figSvg: () => famBandSvg(pick, ref, side), figCss: FAM_CSS,
        subHtml: 'The IDR ends beside ' + esc(pr.name || '') + ' (' + esc(pick || '') + '), ' + esc(FAM_REF.find((x) => x[0] === ref)[1]) + ': excess K/R (blue) and D/E (red) in each depth band, 95% CI by bootstrap over sequence clusters; filled points exclude 0. Click a bar above to show another family.',
        subText: (pr.name || '') + ': K/R and D/E excess by depth from the IDR end, ' + FAM_REF.find((x) => x[0] === ref)[1] + '; 95% CI over sequence clusters.',
        figTsv: () => 'band\tmeasure\texcess_points\tci_low\tci_high\tunits\n' + tableRows(FAM.family_bands).filter((r) => r.interpro === pick && r.reference === ref && r.side === side).map((r) => [r.band, r.measure, 100 * r.excess, 100 * r.ci_low, 100 * r.ci_high, r.units].join('\t')).join('\n') + '\n' };
      body = card(sec.spec) + card(band) + famNotesHtml();
    } else body = card(sec.spec, CTX_ORDER().includes(key) ? '<p class="note">Each modification’s residues: ' + OV.kinds.map((k) => esc(k[1]) + ' ' + k[2].split('').join('/')).join(', ') + '. ADP-ribosylation is tested on serines; ubiquitin and SUMO on lysines.</p>' : '');
    view().innerHTML = selPanel() + nav + '<p class="sec-kicker">' + esc(t[1]) + ' · ' + esc(tileTitle(key)) + '</p><h2 class="sec-q">' + esc(t[2]) + '</h2>' + body + sectionNotes(key) + nav;
    $$('.hbars').forEach((el) => {
      const id = el.closest('.chart-card').id, spec = UI.charts[id];
      if (spec) C.bindBars(el, spec.rows, (r) => { if (sec.families && id.startsWith('ch-fam-')) { UI.famPick = r.key; sectionView('families'); } });
    });
    bindExports(view());
    const on = (id, f) => { const el = $('#' + id); if (el) el.addEventListener('change', () => { f(el.value); sectionView(key); }); };
    on('agm', (v) => { UI.agMeasure = v; }); on('kinctx', (v) => { UI.kinCtx = v; });
    on('famref', (v) => { UI.famRef = v; }); on('famside', (v) => { UI.famSide = v; }); on('famm', (v) => { UI.famMeasure = v; });
    bindSelPanel('overview/' + key);
  }
  function agSummaryHtml() {
    const S = tableRows(OV.agnostic_summary), cp = S.find((r) => r.end === 'domain change points');
    const ends = (k) => (tableRows(OV.agnostic).find((r) => r.end === k && r.measure === 'basic' && r.bin === '1-5') || {}).ends;
    return '<div class="panel note-panel"><p><b>What this shows.</b> Read without drawing a window, the residues beside most folded domains are not a distinct segment: beyond the first five residues (a small basic cap) they look like the rest of their IDR, and protein termini show the same small effects. ' +
      'The exception is the DNA-binding domain of a transcription factor: its ends (' + fmt(ends('domain:dna')) + ') carry a basic extension that fades over about 15 residues, the basic tails known to tune DNA binding. Histone-chaperone folds are too few to read here (' + fmt(ends('domain:chaperone')) + ' ends), with acidic sequence 30–75 residues out. ' +
      (cp ? 'Splitting each domain-touching IDR where its charge changes most finds a split stronger than in shuffles in ' + Math.round(100 * cp.excess_1_5) + '% of ' + fmt(cp.ends) + ' ends (5% expected by chance), at a median of ' + cp.half_depth + ' residues. ' : '') +
      'So a domain flank is an address next to a fold (the dtIDR window) whose content varies by domain class, not a universal element; the tests on this site are tests within that frame. <a href="#overview/families">Domain families</a> reads the same profile for every InterPro entry.</p></div>';
  }
  function famNotesHtml() {
    const both = famRows('idr', 'both'), bothP = famRows('protein', 'both');
    const byK = (a, b) => b.basic_proximal - a.basic_proximal;
    const up = both.filter((r) => SIG(r.basic_proximal_q) && r.basic_proximal > 0).sort(byK), upP = bothP.filter((r) => SIG(r.basic_proximal_q) && r.basic_proximal > 0).sort(byK);
    const acP = bothP.filter((r) => SIG(r.acidic_proximal_q) && r.acidic_proximal > 0), ac = both.filter((r) => SIG(r.acidic_proximal_q) && r.acidic_proximal > 0);
    const only = upP.filter((r) => !up.some((x) => x.interpro === r.interpro)).map((r) => r.name);
    return '<div class="panel note-panel"><p><b>What this shows.</b> Against the rest of their own IDR, ' + fmt(up.length) + ' of ' + fmt(both.length) + ' families have a basic extension at the domain edge: ' + up.map((r) => esc(r.name) + ' (' + pts(r.basic_proximal, 0) + ', ~' + (r.basic_length || '<5') + ' aa)').join(', ') + '. ' +
      'Against the protein’s non-domain residues, ' + fmt(upP.length) + ' do' + (only.length ? ', adding ' + only.map(esc).join(', ') + ': their IDRs are basic throughout rather than at the edge' : '') + '. ' +
      (ac.length + acP.length ? '' : '<b>No family has an acidic extension at the edge</b>, against either reference: acidic flanks are acidic IDRs that happen to touch a fold, not acidity tethered to it. ') +
      'Paralogs are collapsed first: before that, 22 FAM90A copies made a “zinc knuckle” look acidic (+24 points). Small sets (10–30 units) are noisy; read the length with the band chart.</p></div>';
  }
  function sectionNotes(key) {
    const N = {
      acidic_flank: 'Acidic flanks are defined from sequence alone (disorder, acidic charge, glutamate tracts); no modification enters the score. CK2 and glutamylation are the expected marks; read methyl-R and acetyl-K as the basic flank’s marks seen from the other side.',
      basic_flank: 'Basic flanks are the top tenth of MicroDomains’ dtIDR score (disorder, net positive charge, local R/K clustering, RG/GR density); no modification enters it.',
      charge: 'A modification that adds negative charge next to an acidic block lengthens it and raises δ; one inside a basic patch interrupts it and lowers δ. Random placement of as many residues of the same type in the same IDR is the null, so a composition change alone scores zero. Glutamylation is left out: it deepens an existing negative charge, which δ (counted by sign) does not see. Ubiquitin and SUMO attach a protein and are left out.',
      kin: 'Motif matches are predictions: a residue whose neighbours are acidic matches acidophilic kinases whether or not they act on it. The test therefore compares observed phosphosites with other observed phosphosites of the same proteins; it says which kinases could write a context’s sites.',
      crosstalk: 'Methylarginines are enriched at −2 and +2 to +5 from phosphosites, not at the −3 position the basophilic R-x-x-S/T motif needs: shared sequence context (S/G/R-rich segments) more than a kinase-motif link.',
      partners: 'Proteins with acidic or basic flanks are enriched in nuclear and chromatin proteins, which bind each other, so the null keeps degree and nuclear status. Like binds like: basic–basic pairs are over-represented, acidic–basic pairs are not.',
      occupancy: 'Occupancy is the share of a site’s protein copies that carry the modification. Mann–Whitney compares the distributions; read the medians with the counts, as the flanks hold few measured sites.',
      classes: 'Classes come from MicroDomains’ flank table (TFs: Lambert et al. 2018; RNA-binding domains; curated histone chaperones; histone folds; HMG proteins).',
      ctype: 'Two numbers decide the words: net charge per residue (acidic or basic beyond ±0.10, five charges per 50 residues) and, for the rest, the share of residues charged (a quarter, Das & Pappu’s boundary of weakly charged chains). An IDR is the natural unit: the whole disordered segment up to the next domain or the protein’s end, as an acidic tail like HMGB1’s sits beyond a basic linker rather than at the fold. Charge segregation (κ) and the Das–Pappu classes are in the IDR and flank tables.',
    };
    return N[key] ? '<div class="panel note-panel"><p>' + esc(N[key]) + '</p></div>' : '';
  }

  // The ideas the numbers support, each pointing to its tile; numbers from overview.json.
  function ideasHtml() {
    const rt = tableRows(OV.residue_tests), g = (k, c) => rt.find((r) => r.kind === k && r.context === c) || {};
    const or = (k, c) => fx(g(k, c).mh_or) + '×';
    const go = (k, l) => '<a href="#overview/' + k + '">' + esc(l) + ' →</a>';
    const ideas = [
      ['neg', 'Tested: the unit', 'Beside most domains there is no distinct segment', 'Read without a window, the residues next to a fold resemble the rest of their IDR except for a five-residue basic cap; termini show the same. Family by family, basic extensions sit beside DNA-binding folds (homeodomain, forkhead, HMG box, nuclear-receptor zinc finger), for about 5–20 residues; no family has an acidic one.', [go('agnostic', 'Flank profiles'), go('families', 'Domain families')]],
      ['neg', 'Tested: the acidic side', 'Acidity is not tethered to the fold', 'The acidic flanks are acidic IDRs that touch a domain: against either its own IDR or the protein’s other disordered residues, no family is more acidic at the edge. HMGB1 shows the shape: a basic linker against HMG box B, its acidic tail beyond it. Score the whole disordered segment, not a fixed window, for acidic biology.', [go('ctype', 'Charge of IDRs and flanks'), go('families', 'Domain families')]],
      ['data', 'The frame', 'Scored from sequence, tested against modifications', 'Every domain edge of the reviewed proteome is scored from sequence alone, once for basic and once for acidic character; no modification enters either score, so each tile is a test of the scores against data they never saw.', [go('acidic_flank', 'Acidic flanks'), go('basic_flank', 'Basic flanks')]],
      ['new', 'Two chemistries', 'Each charge has its own writers', 'CK2 sites sit in acidic flanks (' + or('ck2', 'acidic_flank') + ') and acidic stretches (' + or('ck2', 'acidic_stretch') + '); methylarginines in RG/RGG clusters (' + or('mR', 'rg_cluster') + ') and basic flanks (' + or('mR', 'basic_flank') + '). The phosphosites of each match different kinases: acidophilic (CK2, GRK, FAM20C) in one, basophilic (PKC, AKT, PKA) in the other. Known cases: CK2 at MYB’s and MAX’s DNA-binding domains, CK2 at HDAC1’s acidic tail.', [go('kin', 'Kinase motifs')]],
      ['new', 'Charge tuning', 'Modifications are placed to change charge blocks', 'Phosphorylation and citrullination at their observed sites sharpen charge blocks more than the same number of random residues would; ADP-ribosylation softens them, as its −2 lands inside basic patches (' + or('adpr', 'basic_patch') + '). Placement, not just composition, is under selection.', [go('charge', 'Charge blocks')]],
      ['new', 'Arginine', 'PADs and PRMTs pick different arginines', 'Citrullines avoid RG/RGG clusters (' + or('cit', 'rg_cluster') + ') where methylarginines crowd (' + or('mR', 'rg_cluster') + '); both favour IDRs. The two arginine marks compete for few residues.', [go('rg_cluster', 'RG/RGG clusters')]],
      ['neg', 'Tested: not a rule', 'Acidic and basic flanks do not pair across partners', 'Opposite-charge partner pairs are not over-represented after degree and nuclear status are held; same-charge pairs are. Charge pairing is real in named systems (chaperone and histone), not who binds whom across the proteome.', [go('partners', 'Partners')]],
      ['hyp', 'Hypothesis', 'A charged flank tunes avidity', 'The folded domain chooses the site; the adjacent residues make a second, low-specificity contact with the same ligand surface. A modification placed to sharpen or soften a charge block retunes that contact, so it should change residence time rather than site selection.', [go('charge', 'Charge blocks')]],
      ['open', 'Open', 'What is not settled', 'Whether CK2 writes the acidic IDRs it matches or only marks them; why phosphosites in acidic flanks are not more occupied; a proteome-wide glutamylation site table (none exists; NAP1L1 and NAP1L4 are glutamylated inside their acidic tails); FINCHES predictions of IDR–IDR attraction for the partner test.', []],
    ];
    return '<section class="ideas-sec"><h2>Ideas</h2><div class="ideas">' + ideas.map(([tag, tl, h, p, links]) => '<article class="idea"><span class="tag ' + tag + '">' + esc(tl) + '</span><h3>' + esc(h) + '</h3><p>' + p + '</p>' + (links.length ? '<div class="foot">' + links.join('') + '</div>' : '') + '</article>').join('') + '</div></section>';
  }

  // ---------------------------------------------------------------- Definitions
  async function definitionsView() {
    const ctx = OV.contexts.map(([k, l, d]) => '<dt id="def-ctx-' + k + '">' + esc(l) + '</dt><dd>' + esc(d) + '.</dd>').join('');
    const kinds = OV.kinds.map(([k, l, res, mod]) => '<tr><td>' + esc(l) + '</td><td>' + res.split('').join(', ') + '</td><td>' + esc(mod || '—') + '</td><td class="num">' + fmt(B.counts.by_kind[k] || 0) + '</td></tr>').join('');
    const S = OV.settings;
    view().innerHTML = '<div class="about defs"><h2>Definitions</h2>' +
      '<div class="panel"><h3 id="def-unit">Domain flanks</h3><dl>' +
      '<dt>Domain flank</dt><dd>The residues beside one edge of a folded domain (InterPro domain or repeat), up to 50 on that side or fewer where the next domain begins: the window a domain-tethered IDR (dtIDR) would occupy. ' + fmt(B.counts.flanks) + ' flanks on ' + fmt(B.counts.proteins_with_flanks) + ' proteins. N-side: before the domain; C-side: after it.</dd>' +
      '<dt>Basic score</dt><dd>MicroDomains’ dtIDR score (<code>dtidr_score_z</code>): disorder, net positive charge, local R/K clustering and RG/GR density, each as a percentile, averaged and standardised. Shown as a z-score.</dd>' +
      '<dt>Acidic score</dt><dd>The same construction with the charge inverted: disorder, acidic charge fraction and density of glutamates in acidic runs (<code>acidic_dtidr_score</code>).</dd>' +
      '<dt id="def-top">Acidic flank, basic flank</dt><dd>A flank in the top tenth of all flanks (' + (100 * S.flank_top).toFixed(0) + 'th percentile) for that score: ' + fmt(B.counts.acidic_flanks) + ' acidic and ' + fmt(B.counts.basic_flanks) + ' basic flanks. A flank can be both. Earlier versions called these the acidic and basic poles.</dd>' +
      '<dt>IDR</dt><dd>A MobiDB disordered region (overlapping fragments merged) of ' + S.min_idr + ' or more residues: ' + fmt(B.counts.idrs) + '. The whole disordered segment, up to the next domain or the protein’s end.</dd></dl></div>' +
      '<div class="panel"><h3 id="def-ctype">Charge, in words</h3><p>Every IDR and flank is called one of four things, from two numbers anyone can check on its sequence (idrkit.charge_type):</p><dl>' +
      [3, 4, 2, 1].map((t) => '<dt>' + ctChip(t) + '</dt><dd>' + esc(CT_DEF[t]) + '.</dd>').join('') +
      '<dt>Net charge</dt><dd>K and R count +1, D and E −1, histidine 0; summed over the segment, or per residue (NCPR). The cut of ±0.10 per residue is five charges in excess per 50 residues.</dd>' +
      '<dt>Charged share</dt><dd>D, E, K and R over all residues (FCR). The quarter is Das & Pappu’s boundary between weakly and strongly charged chains (PNAS 110:13392, 2013).</dd>' +
      '<dt>On the plots</dt><dd>% acidic residues across, % basic residues up: below the red line acidic, above the blue line basic, between them mixed or, near the corner, few charges.</dd></dl></div>' +
      '<div class="panel"><h3 id="def-cider">Charge patterning (CIDER)</h3><p>Every IDR and every flank is also described by localCIDER 0.1.21 (Holehouse et al., Biophys J 112:16, 2017), the parameters of Das & Pappu (PNAS 110:13392, 2013) and Martin et al. (2016), without any score of ours:</p><dl>' +
      '<dt>FCR, NCPR, f+, f−</dt><dd>Fractions of charged residues, net charge per residue, fractions of positive (K, R) and negative (D, E) residues.</dd>' +
      '<dt>Diagram-of-states class</dt><dd>Das & Pappu’s five classes predict a disordered chain’s shape from its fractions of positive and negative residues; each is named here in plain words with the field’s term beside it:<ul class="dp-list">' + [1, 2, 3, 4, 5].map((r) => '<li>' + dpChip(r) + ' ' + esc(REGION[r]) + ': ' + esc(REGION_RULE[r]) + '; predicted shape <b>' + esc(REGION_SHAPE[r]) + '</b></li>').join('') + '</ul></dd>' +
      '<dt>κ (kappa)</dt><dd>How segregated the positive and negative charges are: δ divided by the largest δ the same residues can reach; 0 well mixed, 1 fully segregated. Left empty for fewer than two charges of each sign.</dd>' +
      '<dt>Ω (Omega)</dt><dd>κ’s extension to charged residues and prolines together.</dd>' +
      '<dt>Hydropathy</dt><dd>Mean Kyte–Doolittle hydropathy normalised to 0–1.</dd></dl><p class="note">idrkit’s δ is checked against localCIDER’s for every IDR and flank when the parameters are computed.</p></div>' +
      '<div class="panel"><h3 id="def-agnostic">Reading a flank without a window</h3><p>Every IDR is read residue by residue from each of its ends. At each depth, a measure (K/R fraction, D/E fraction, net charge, the |net charge| of the next 10 residues, or the share of residues of a type that carry a modification) is compared with the same IDR’s residues in random order (20 shuffles), so length and composition are kept and only position is removed. Ends are classed by what lies beyond them: a folded domain (InterPro domain or repeat within 10 residues; split by the domain’s class), the protein’s terminus, or ordered sequence that is not a domain. 95% intervals: bootstrap over proteins. A per-end change point (the split with the largest proximal–distal difference in charge concentration, against shuffles) gives the length scale without choosing one.</p></div>' +
      '<div class="panel"><h3 id="def-families">Domain families</h3><p>The same reading for every InterPro domain and repeat separately: each IDR end that touches the entry (within 10 residues), its K/R and D/E fractions and net charge in depth bands 1–5, 6–10, 11–15, 16–20, 21–30 and 31–50, minus a reference: the same IDR (is the edge a distinct segment?) or the protein’s residues outside domains (is charge tethered to the fold at all?). A protein’s ends are averaged, then proteins whose ends share an identical 12-residue stretch within 30 residues of the end (paralogs) are averaged into one sequence cluster, the unit tested. One-sample t-test of residues 1–10 and 11–30 against 0, BH within a reference; an extension’s length is the deepest band, contiguous from the edge, whose 95% bootstrap interval excludes 0. Entries with ends on at least 10 clusters.</p></div>' +
      '<div class="panel"><h3 id="def-contexts">Where a residue sits</h3><p>Each test asks whether a modification’s residues sit in a context more often than the unmodified residues of the same type in the same proteins. A residue can be in several contexts; each is tested on its own, inside against outside.</p><dl>' + ctx + '</dl></div>' +
      '<div class="panel"><h3 id="def-kinds">Modifications</h3><p>Sites as the sources report them, placed on the UniProt 2026_03 reviewed sequences by residue and, where the source gives peptides, by peptide. The charge model is the one the charge tests use (idrkit.ptm).</p>' +
      '<div class="table-wrapper" tabindex="0" role="region" aria-label="Modification kinds"><table class="data"><thead><tr><th>Modification</th><th>Residues tested</th><th>Charge model</th><th class="num">Sites</th></tr></thead><tbody>' + kinds + '</tbody></table></div>' +
      '<p class="note">Charge of a modified residue at pH 7.4: phospho −2; acetyl-K and neutral acyls 0; succinyl, malonyl and glutaryl-K −1; citrulline 0; ADP-ribosyl-S −2 (on D or E: −2 in place of −1); one added glutamate −2 in place of −1 (a floor); methylation unchanged. Ubiquitin and SUMO are left out of charge.</p></div>' +
      '<div class="panel"><h3 id="def-charge">Charge patterning</h3><dl>' +
      '<dt>δ (delta)</dt><dd>Das & Pappu’s charge patterning: the mean over all blobs of 5 and of 6 residues of (σ − σ<sub>blob</sub>)², σ = NCPR²/FCR, as localCIDER computes it (idrkit.delta). Higher: charges sit in local blocks.</dd>' +
      '<dt>δ z</dt><dd>δ against 40 shuffles of the same residues (length and composition kept). Above 0: blockier than the composition gives. An IDR is <b>blocky</b> at δ z ≥ ' + S.blocky_z + '. κ normalises δ by its maximum; δ z by its spread under shuffles.</dd>' +
      '<dt>SCD z</dt><dd>Sequence charge decoration (Sawle & Ghosh 2015) against the same shuffles: long-range segregation; below 0, more segregated.</dd>' +
      '<dt>Sharper or softer than random</dt><dd>For each IDR and modification, δ with the observed sites modified against δ with as many residues of the same type in the same IDR modified at random (20 draws). The share of IDRs where the observed change is larger is tested with a sign test.</dd></dl></div>' +
      '<div class="panel"><h3 id="def-kinases">Kinase motifs</h3><p>Every serine and threonine in an IDR or flank, and every observed phosphosite, scored with the Ser/Thr kinase matrices of Johnson et al. 2023 (Nature 613:759) through kinase-library 1.8.0: for each kinase, the percentile of the residue’s score among the reference human phosphoproteome. A <b>motif match</b> is a percentile of ' + OV.kinase_match + ' or more. Representatives name each motif class: ' +
      OV.kinase_rows.map((r) => esc(r[1]) + ' (' + esc(r[2]) + ')').join(', ') + '. Strong unobserved CK2 motifs (CK2 percentile ≥ 99) are drawn on the protein maps: ' + fmt(B.counts.dark_ck2_motifs) + '.</p></div>' +
      '<div class="panel"><h3 id="def-occupancy">Occupancy</h3><p>Phospho-S/T occupancy: the share of a protein’s copies carrying the phosphate, measured in 33 cancer cell lines (Koenig et al. 2026), median over the lines, unstimulated. Acetyl-K stoichiometry: the same for acetylation, measured in HeLa (Hansen et al. 2019).</p></div>' +
      '<div class="panel"><h3 id="def-tests">Tests</h3><dl>' +
      '<dt>Within protein</dt><dd>Mantel–Haenszel common odds ratio, each protein a stratum (its sites and its unmodified residues of the same type), with the Robins–Breslow–Greenland 95% CI and the continuity-corrected CMH test (PTM_Pipeline <code>common.mantel_haenszel_ci</code>; in the browser, the Methylarginine Site Explorer’s <code>mhTest</code>). It asks whether, within a protein, the modified residues lie in a context more often than its other residues, which a pooled test cannot separate from which proteins are modified. A test needs ' + S.min_in_context + ' residues in the context and ' + S.min_strata + ' informative proteins.</dd>' +
      '<dt>Pooled</dt><dd>Fisher’s exact test with the Haldane–Anscombe odds ratio over the same counts (shown on hover); the unit for the Das–Pappu classes, where IDRs are counted.</dd>' +
      '<dt>Protein classes</dt><dd>One-sided hypergeometric test, fold enrichment (k/n ÷ M/N) against all reviewed human proteins.</dd>' +
      '<dt>Distributions</dt><dd>Mann–Whitney U, normal approximation with tie and continuity corrections.</dd>' +
      '<dt>Multiple testing</dt><dd>Benjamini–Hochberg across the rows of each chart; filled marks have q ≤ ' + PADJ + '. Effects within ' + C.minEffect() + '× either way are drawn grey (about the same) whatever their q.</dd>' +
      '<dt>Your selection</dt><dd>The same tests run again in the browser on the selected proteins and evidence, from the same per-protein counts; sites the selection does not keep are left out of both sides. The tick on each bar is the whole dataset. A selection’s q values are exploratory: they are not adjusted for the other selections tried.</dd></dl></div>' +
      '<div class="panel"><h3>How to read the figures</h3><p>Tiles and charts are the Methylarginine Site Explorer’s, drawn by its own code: a dot (or bar) at the odds ratio on a log scale, filled red where more often, slate where less, grey where about the same, hollow where not significant; a line through a dot is its 95% CI. A chart’s title names what it plots; the question it answers heads its section. Charts export as SVG, PNG and TSV; protein maps as SVG and PNG of the view shown.</p></div></div>';
  }

  // ---------------------------------------------------------------- Proteins, Domain flanks, IDRs, Sites
  const inProtSel = () => (R.scoped && R.sel && R.sel.prots ? R.sel.prots : null);
  const KCOLS = (D, on) => OV.kinds.map(([k, l, res]) => ({ k, label: l, num: true, opt: true, on: on.includes(k), menu: l + ' sites', get: (i) => D[k][i], tip: l + ' sites (' + res.split('').join('/') + ')' }));
  const siteTotal = (D, i) => OV.kinds.reduce((s, [k]) => s + (k === 'ck2' ? 0 : D[k][i] || 0), 0);
  const geneCell = (p, extra) => '<span class="gene">' + esc(geneOf(p)) + '</span><span class="acc">' + esc(PR.acc[p]) + (extra ? ' · ' + extra : '') + '</span>';
  const TOPB = (bits) => (bits & 2 ? ' <span class="top b" title="basic flank: top tenth of all flanks by the basic score">B</span>' : '') + (bits & 1 ? ' <span class="top a" title="acidic flank: top tenth of all flanks by the acidic score">A</span>' : '');
  const qMatch = (q) => { const Q = q.trim().toUpperCase(); return Q ? (p) => (PR.gene[p] || '').toUpperCase().startsWith(Q) || PR.acc[p].startsWith(Q) : null; };
  const pctTop = (v) => (v == null ? '' : v >= 990 ? 'top 1%' : v >= 900 ? 'top ' + Math.ceil((1000 - v) / 10) + '%' : Math.round(v / 10) + 'th pct');
  // the filters above a table: a query box, selects and checks, in the URL; a change re-routes
  function filterForm(id, st, items) {
    return '<form class="filters" id="' + id + '">' + items.map((it) => {
      if (it.type === 'text') return '<label>' + esc(it.label) + '<input type="' + (it.num ? 'number' : 'search') + '" name="' + it.k + '" value="' + esc(st[it.k] || '') + '" placeholder="' + esc(it.ph || '') + '"' + (it.num ? ' min="0" step="10"' : '') + '></label>';
      if (it.type === 'select') return '<label>' + esc(it.label) + '<select name="' + it.k + '"><option value="">' + esc(it.any || 'any') + '</option>' + it.opts.map(([v, l]) => '<option value="' + v + '"' + (String(st[it.k]) === String(v) ? ' selected' : '') + '>' + esc(l) + '</option>').join('') + '</select></label>';
      return '<label class="check"><input type="checkbox" name="' + it.k + '"' + (st[it.k] ? ' checked' : '') + '>' + esc(it.label) + '</label>';
    }).join('') + '</form>';
  }
  const readFilters = (params, items) => Object.fromEntries(items.map((it) => [it.k, it.type === 'check' ? params.get(it.k) === '1' : params.get(it.k) || '']));
  const CT_OPTS = [[3, 'Acidic'], [4, 'Basic'], [2, 'Mixed charge'], [1, 'Few charges']];
  // A table view: selection, filters, an optional plot that selects, the table. o.plot: { axes: {key: {label, get,
  // fmt, log}}, ax, ay (defaults), colours: {key: {label, cats, cat}}, col, guides(ax, ay) }
  async function tableView(o) {
    const { path, params, items, filter, n, cols, sortKey, sortDir, plot } = o;
    const st = readFilters(params, items);
    st.tbl = tblState(params, sortKey, sortDir);
    const P = plot ? { ax: plot.axes[params.get('ax')] ? params.get('ax') : plot.ax, ay: plot.axes[params.get('ay')] ? params.get('ay') : plot.ay, col: plot.colours[params.get('col')] ? params.get('col') : plot.col,
      box: (params.get('box') || '').split(',').filter(Boolean).map(Number) } : null;
    if (P && (P.box.length !== 4 || P.box.some((v) => !Number.isFinite(v)))) P.box = null;
    const pass = [];
    for (let i = 0; i < n; i++) if (filter(i, st)) pass.push(i);
    const AX = P && plot.axes[P.ax], AY = P && plot.axes[P.ay];
    const inBox = P && P.box ? (i) => { const x = AX.get(i), y = AY.get(i); return x != null && y != null && x >= P.box[0] && x <= P.box[1] && y >= P.box[2] && y <= P.box[3]; } : null;
    const rows = inBox ? pass.filter(inBox) : pass;
    const save = () => {
      const p = selParams();
      for (const it of items) { const v = st[it.k]; if (it.type === 'check' ? v : v !== '' && v != null) p.set(it.k, it.type === 'check' ? '1' : v); }
      tblParams(p, st.tbl, sortKey, sortDir);
      if (P) { if (P.ax !== plot.ax) p.set('ax', P.ax); if (P.ay !== plot.ay) p.set('ay', P.ay); if (P.col !== plot.col) p.set('col', P.col); if (P.box) p.set('box', P.box.join(',')); }
      setParams(path, p);
    };
    const opt = (id, label, map, cur) => '<label>' + esc(label) + ' <select id="' + id + '">' + Object.entries(map).map(([k, v]) => '<option value="' + k + '"' + (k === cur ? ' selected' : '') + '>' + esc(v.label) + '</option>').join('') + '</select></label>';
    view().innerHTML = selPanel() + C.resultsHead(esc(o.title), '<b>' + fmt(rows.length) + '</b> of ' + fmt(n) + ' ' + esc(o.unit) + (P && P.box ? ' · selected on the plot' : '') + (R.scoped ? ' · ' + esc(selectionText()) : ''), '') +
      '<p class="lede">' + o.lede + '</p><div class="panel">' + filterForm('ff', st, items) + '</div>' +
      (plot ? '<div class="card plot-card"><div class="card-head"><div><h3>' + esc(AY.label) + ' against ' + esc(AX.label.replace(/^[A-Z](?=[a-z])/, (c) => c.toLowerCase())) + '</h3><p class="sub">' + esc(plot.sub) + '</p></div></div>' +
        '<div class="chart-tools">' + opt('pax', 'Across', plot.axes, P.ax) + opt('pay', 'Up', plot.axes, P.ay) + opt('pcol', 'Colour by', plot.colours, P.col) + '</div><div id="tplot"></div></div>' : '') + '<div id="tt"></div>';
    const form = $('#ff');
    form.addEventListener('change', (e) => { const t = e.target; st[t.name] = t.type === 'checkbox' ? t.checked : t.value; st.tbl.n = PAGE_MORE; save(); route(); });
    form.addEventListener('submit', (e) => e.preventDefault());
    const draw = () => tbl($('#tt'), { id: path, cols, rows, st: st.tbl, caption: o.caption, note: o.note, tsv: o.tsv, onOpen: o.open, onState: () => { save(); draw(); } });
    draw();
    if (plot) {
      const C2 = plot.colours[P.col];
      scatter($('#tplot'), { n, x: AX.get, y: AY.get, xLabel: AX.label, yLabel: AY.label, xFmt: AX.fmt, yFmt: AY.fmt, xLog: !!AX.log, rows: pass, inBox, cats: C2.cats, cat: C2.cat, tip: plot.tip,
        onOpen: o.open, box: P.box, onBox: (b) => { P.box = b; st.tbl.n = PAGE_MORE; save(); route(); }, guides: plot.guides && plot.guides(P.ax, P.ay, P.col) });
      for (const [id, k] of [['pax', 'ax'], ['pay', 'ay'], ['pcol', 'col']]) $('#' + id).addEventListener('change', (e) => { P[k] = e.target.value; if (k !== 'col') P.box = null; save(); route(); });
    }
    bindSelPanel(path);
  }
  // the charge map's guides: net charge of 0.10 per residue either way, and a quarter of residues charged
  // with Das-Pappu colours, its boundaries instead: a quarter and a third of residues charged, and a third of one sign
  const chargeGuides = (ax, ay, col) => ((ax === 'fm' && ay === 'fp') ? (g, X, Y, W, H, M) => {
    const ln = (pts, colour) => { g.strokeStyle = colour; g.setLineDash([5, 4]); g.beginPath(); pts.forEach(([a, b], k) => (k ? g.lineTo(X(a), Y(b)) : g.moveTo(X(a), Y(b)))); g.stroke(); g.setLineDash([]); };
    const c = (v) => C.cssVar(v);
    if (col === 'dp') {
      ln([[0, 25], [25, 0]], c('--muted')); ln([[0, 35], [35, 0]], c('--muted')); ln([[0, 35], [65, 100]], c('--dp5')); ln([[35, 0], [100, 65]], c('--dp4'));
      g.font = '600 12px ' + (c('--font-body') || 'sans-serif'); g.textAlign = 'left';
      g.lineWidth = 3; g.strokeStyle = c('--surface'); g.lineJoin = 'round';
      for (const [r, a, b] of [[1, 2, 3], [2, 22, 8], [3, 28, 24], [4, 50, 6], [5, 2, 44]]) { const t = 'R' + r + ' ' + REGION_SHORT[r].slice(3); g.strokeText(t, X(a), Y(b)); g.fillStyle = c('--dp' + r); g.fillText(t, X(a), Y(b)); }
      g.lineWidth = 1;
      return;
    }
    ln([[10, 0], [100, 90]], c('--acidic')); ln([[0, 10], [90, 100]], c('--basic')); ln([[0, 25], [25, 0]], c('--muted'));
    g.font = '600 12px ' + (c('--font-body') || 'sans-serif'); g.textAlign = 'left';
    const lab = (t, a, b, colour) => { g.lineWidth = 3; g.strokeStyle = c('--surface'); g.lineJoin = 'round'; g.strokeText(t, X(a), Y(b)); g.lineWidth = 1; g.fillStyle = colour; g.fillText(t, X(a), Y(b)); };
    lab('acidic', 30, 4, c('--acidic')); lab('basic', 2, 32, c('--basic')); lab('mixed', 22, 20, c('--hybrid')); lab('few charges', 2, 3, c('--muted'));
  } : null);
  const pctAxis = (label, get) => ({ label, get, fmt: (v) => v + '%' });

  async function proteinsView(params) {
    view().innerHTML = '<p class="loading">Loading 20k proteins…</p>';
    await Promise.all([loadProteins(), computeSelection()]);
    const ps = inProtSel(), cls = (p) => OV.classes.filter((c, b) => PR.cls[p] & (1 << b)).map((c) => c[1]);
    const total = (p) => siteTotal(PR, p);
    const mix = (p) => { const t = total(p); if (!t) return '<span class="muted">—</span>'; return '<span class="mix" data-tip="' + esc(OV.kinds.filter(([k]) => k !== 'ck2' && PR[k][p]).map(([k, l]) => fmt(PR[k][p]) + ' ' + l).join('<br>')) + '">' + OV.kinds.filter(([k]) => k !== 'ck2' && PR[k][p]).map(([k]) => '<i style="flex-grow:' + PR[k][p] + ';background:var(--' + KCOL[k] + ')"></i>').join('') + '</span><b class="mix-n">' + fmt(t) + '</b>'; };
    const ctCell = (p) => [3, 4, 2].filter((t) => PR.ct[p][t - 1]).map((t) => PR.ct[p][t - 1] + ' ' + ctChip(t)).join(' ') || '<span class="muted">' + (PR.nidr[p] ? 'few charges' : 'no IDR') + '</span>';
    const items = [{ type: 'text', k: 'q', label: 'Gene or accession', ph: 'e.g. NPM1' }, { type: 'select', k: 'ct', label: 'Has an IDR that is', opts: CT_OPTS },
      { type: 'check', k: 'top', label: 'an acidic or basic flank' }, { type: 'check', k: 'mod', label: 'any modification site' }];
    await tableView({ path: 'proteins', params, title: 'Proteins', unit: 'reviewed human proteins', items, n: PR.acc.length, sortKey: 'sites', sortDir: -1, caption: 'Proteins', tsv: 'idr_explorer_proteins.tsv',
      lede: 'Every reviewed human protein: how much of it is disordered, how its disordered regions are charged, its best domain flank of each charge and its modification sites. Click a protein to open it.' + docLink('unit'),
      filter: (p, st) => (!ps || ps.has(p)) && (!st.q || qMatch(st.q)(p)) && (!st.ct || PR.ct[p][+st.ct - 1] > 0) && (!st.top || PR.has[p] & 3) && (!st.mod || total(p) > 0),
      open: (p) => { location.hash = PR.acc[p]; },
      cols: [
        { k: 'gene', label: 'Protein', get: (p) => geneOf(p), td: (p) => geneCell(p) + '<div class="pname">' + esc(PR.name[p]) + '</div>' },
        { k: 'len', label: 'Length', num: true, get: (p) => PR.len[p] },
        { k: 'dis', label: 'Disordered', num: true, get: (p) => PR.dis[p], td: (p) => PR.dis[p] + '%', tip: 'Share of residues in MobiDB disordered regions' },
        { k: 'nidr', label: 'IDRs', num: true, get: (p) => PR.nidr[p], tip: 'Disordered regions of ' + OV.settings.min_idr + ' or more residues' },
        { k: 'ct', label: 'Charged IDRs', get: (p) => PR.ct[p][2] * 1000 + PR.ct[p][3] * 10 + PR.ct[p][1], td: ctCell, tip: 'Its IDRs that are acidic, basic or of mixed charge (net charge per residue and fraction charged)' },
        { k: 'dp', label: 'Das–Pappu classes', opt: true, get: (p) => PR.dp[p] || null, td: (p) => dpChips(PR.dp[p]), menu: 'Das–Pappu classes of its IDRs' },
        { k: 'za', label: 'Best acidic flank', num: true, get: (p) => z(PR.za[p]), td: (p) => zChip(PR.za[p], 'acidic') + TOPB(PR.has[p] & 1), tip: 'Highest acidic score (z) of its domain flanks; A: in the top tenth of all flanks' },
        { k: 'zb', label: 'Best basic flank', num: true, get: (p) => z(PR.zb[p]), td: (p) => zChip(PR.zb[p], 'basic') + TOPB(PR.has[p] & 2), tip: 'Highest basic score (z) of its domain flanks; B: in the top tenth of all flanks' },
        { k: 'sites', label: 'Modification sites', num: true, get: total, td: mix, tip: 'All modification sites, coloured by kind (hover for the counts)' },
        { k: 'cls', label: 'Class', get: (p) => cls(p).join(', ') || null, td: (p) => cls(p).map((c) => '<span class="tag">' + esc(c) + '</span>').join(' ') },
        { k: 'nas', label: 'Acidic stretches', num: true, opt: true, get: (p) => PR.nas[p] }, { k: 'nbp', label: 'Basic patches', num: true, opt: true, get: (p) => PR.nbp[p] },
        { k: 'nrg', label: 'RG/RGG clusters', num: true, opt: true, get: (p) => PR.nrg[p] },
      ].concat(KCOLS(PR, [])),
    });
  }

  async function flanksView(params) {
    view().innerHTML = '<p class="loading">Loading 121k domain flanks…</p>';
    await Promise.all([loadProteins(), loadFlanks(), computeSelection()]);
    const ps = inProtSel();
    const fm = (i) => (FL.fcr[i] == null ? null : (FL.fcr[i] - FL.ncpr[i]) / 20), fp = (i) => (FL.fcr[i] == null ? null : (FL.fcr[i] + FL.ncpr[i]) / 20);
    const items = [{ type: 'text', k: 'q', label: 'Gene or accession', ph: 'e.g. NCL' },
      { type: 'select', k: 'ct', label: 'Charge', opts: CT_OPTS },
      { type: 'select', k: 'top', label: 'Top tenth by score', opts: [['a', 'acidic flank'], ['b', 'basic flank'], ['ab', 'both']] },
      { type: 'select', k: 'side', label: 'Side of the domain', any: 'either', opts: [['N', 'before it (N-side)'], ['C', 'after it (C-side)']] },
      { type: 'check', k: 'tf', label: 'TF' }, { type: 'check', k: 'rbd', label: 'RNA-binding domain' }, { type: 'check', k: 'chap', label: 'histone chaperone or HMG' },
      { type: 'check', k: 'mr', label: 'methyl-R inside' }, { type: 'check', k: 'ck2', label: 'CK2 site inside' }, { type: 'check', k: 'hist', label: 'include histone-like' }];
    const dn = (i) => (FL.domain_names[FL.dn[i]] || '').split(';')[0];
    const T = (i) => FL.top[i];
    await tableView({ path: 'flanks', params, title: 'Domain flanks', unit: 'flanks', items, n: FL.p.length, sortKey: 'za', sortDir: -1, caption: 'Domain flanks', tsv: 'idr_explorer_domain_flanks.tsv',
      lede: 'The up to 50 residues beside every folded-domain edge (where a domain-tethered IDR would sit), with how charged they are and their two scores. The top tenth by each score are the <b>acidic</b> and <b>basic flanks</b>. Pick flanks on the plot or with the filters; click one to open its protein there.' + docLink('unit'),
      filter: (i, st) => {
        const f = FL.fl[i];
        if (!st.hist && (f & 4)) return false;
        if (ps && !ps.has(FL.p[i])) return false;
        if (st.q && !qMatch(st.q)(FL.p[i])) return false;
        if (st.ct && FL.ct[i] !== +st.ct) return false;
        if (st.side && FL.side[i] !== (st.side === 'N' ? 0 : 1)) return false;
        if (st.top === 'a' && !(T(i) & 1)) return false;
        if (st.top === 'b' && !(T(i) & 2)) return false;
        if (st.top === 'ab' && T(i) !== 3) return false;
        if (st.tf && !(f & 1)) return false;
        if (st.rbd && !(f & 2)) return false;
        if (st.chap && !(f & 24)) return false;
        if (st.mr && !FL.mR[i]) return false;
        if (st.ck2 && !FL.ck2[i]) return false;
        return true;
      },
      open: (i) => { location.hash = PR.acc[FL.p[i]] + '?f=' + FL.fs[i] + '-' + FL.fe[i]; },
      plot: {
        sub: 'Each flank a point: grey all of them, coloured those your filters keep. Hover for the flank, click to open it, drag to select flanks for the table.',
        ax: 'fm', ay: 'fp', col: 'top',
        axes: { fm: pctAxis('% acidic residues (D, E)', fm), fp: pctAxis('% basic residues (K, R)', fp), ncpr: { label: 'Net charge per residue', get: (i) => (FL.ncpr[i] == null ? null : FL.ncpr[i] / 1000) },
          pa: { label: 'Acidic score, percentile', get: (i) => (FL.pa[i] == null ? null : FL.pa[i] / 10) }, pb: { label: 'Basic score, percentile', get: (i) => (FL.pb[i] == null ? null : FL.pb[i] / 10) },
          kap: { label: 'Charge segregation (κ)', get: (i) => (FL.kap[i] == null ? null : FL.kap[i] / 1000) }, sites: { label: 'Modification sites in it', get: (i) => siteTotal(FL, i) } },
        colours: { top: { label: 'acidic or basic flank', cats: [['Acidic flank', '--acidic'], ['Basic flank', '--basic'], ['Both', '--hybrid'], ['Neither', '--ink-2']], cat: (i) => (T(i) === 3 ? 2 : T(i) & 1 ? 0 : T(i) & 2 ? 1 : 3) },
          ct: { label: 'charge', cats: [['Acidic', '--acidic'], ['Basic', '--basic'], ['Mixed charge', '--hybrid'], ['Few charges', '--ink-2']], cat: (i) => [3, 3, 2, 0, 1][FL.ct[i]] },
          dp: { label: 'Das–Pappu class', cats: [1, 2, 3, 4, 5].map((r) => [REGION_SHORT[r], '--dp' + r]), cat: (i) => (FL.reg[i] || 1) - 1 } },
        guides: chargeGuides,
        tip: (i) => '<b>' + esc(geneOf(FL.p[i])) + ' ' + FL.fs[i] + '–' + FL.fe[i] + '</b> (' + (FL.side[i] ? 'after' : 'before') + ' ' + esc(dn(i).slice(0, 40)) + ')<br>' + CT[FL.ct[i]] + ': D/E ' + Math.round(fm(i)) + '%, K/R ' + Math.round(fp(i)) + '%, net ' + fsign(FL.ncpr[i] / 1000, 2) + ' per residue' +
          '<br>acidic score ' + pctTop(FL.pa[i]) + ', basic ' + pctTop(FL.pb[i]) + (siteTotal(FL, i) ? '<br>' + fmt(siteTotal(FL, i)) + ' modification sites' : ''),
      },
      cols: [
        { k: 'gene', label: 'Protein', get: (i) => geneOf(FL.p[i]), td: (i) => geneCell(FL.p[i]) },
        { k: 'fs', label: 'Flank', tdCls: 'nw', get: (i) => FL.fs[i], td: (i) => FL.fs[i] + '–' + FL.fe[i] + ' <span class="muted">' + (FL.side[i] ? 'after' : 'before') + '</span>', tip: 'Residues of the flank, before (N-side) or after (C-side) the domain' },
        { k: 'dn', label: 'Domain', get: (i) => dn(i) || null, td: (i) => '<span class="dn">' + esc(dn(i).slice(0, 44)) + '</span> <span class="muted">' + FL.ds[i] + '–' + FL.de[i] + '</span>' },
        { k: 'ncpr', label: 'Charge', num: true, get: (i) => (FL.ncpr[i] == null ? null : FL.ncpr[i] / 1000), td: (i) => chargeCell(fm(i) / 100, fp(i) / 100, FL.ncpr[i] / 1000, FL.ct[i]), tip: 'D/E (red) and K/R (blue) share of residues, and net charge per residue; sorts by net charge' },
        { k: 'dp', label: 'Das–Pappu class', opt: true, get: (i) => FL.reg[i] || null, td: (i) => dpChip(FL.reg[i]), menu: 'Das–Pappu class (predicted shape)' },
        { k: 'za', label: 'Acidic score', num: true, get: (i) => z(FL.za[i]), td: (i) => '<span class="pt">' + pctTop(FL.pa[i]) + '</span>' + TOPB(T(i) & 1), tip: 'Acidic score (MicroDomains: disorder, acidic charge, glutamate tracts) as a percentile of all flanks; A: top tenth' },
        { k: 'zb', label: 'Basic score', num: true, get: (i) => z(FL.zb[i]), td: (i) => '<span class="pt">' + pctTop(FL.pb[i]) + '</span>' + TOPB(T(i) & 2), tip: 'Basic score (MicroDomains dtIDR score) as a percentile of all flanks; B: top tenth' },
        { k: 'kap', label: 'Charge segregation (κ)', num: true, opt: true, get: (i) => (FL.kap[i] == null ? null : FL.kap[i] / 1000), td: (i) => f3(FL.kap[i]), tip: 'κ (localCIDER): 0 positive and negative residues well mixed, 1 fully segregated; needs two charges of each sign' },
        { k: 'chem', label: 'Chemistry', opt: true, get: (i) => CHEM_SHORT[FL.chem[i]] || null },
        { k: 'lar', label: 'Longest D/E run', num: true, opt: true, get: (i) => FL.lar[i] },
        { k: 'flags', label: 'Class', get: (i) => flagText(FL.fl[i]) || null },
        { k: 'sites', label: 'Sites', num: true, get: (i) => siteTotal(FL, i), tip: 'Modification sites inside the flank (Columns: by kind)' },
      ].concat(KCOLS(FL, ['ph', 'mR'])),
    });
  }

  async function idrsView(params) {
    view().innerHTML = '<p class="loading">Loading 33k IDRs…</p>';
    await Promise.all([loadProteins(), loadIdrs(), computeSelection()]);
    const ps = inProtSel();
    const items = [{ type: 'text', k: 'q', label: 'Gene or accession', ph: 'e.g. NPM1' }, { type: 'text', k: 'minlen', label: 'Minimum length', num: true, ph: String(OV.settings.min_idr) },
      { type: 'select', k: 'ct', label: 'Charge', opts: CT_OPTS },
      { type: 'select', k: 'term', label: 'Position', opts: [['N', 'N-terminal'], ['C', 'C-terminal'], ['I', 'internal']] },
      { type: 'check', k: 'fl', label: 'touches a domain' }, { type: 'check', k: 'as', label: 'acidic stretch' }, { type: 'check', k: 'bp', label: 'basic patch' }, { type: 'check', k: 'rg', label: 'RG/RGG cluster' },
      { type: 'check', k: 'blocky', label: 'charge blocks' }, { type: 'check', k: 'mr', label: 'methyl-R' }, { type: 'check', k: 'ck2', label: 'CK2 site' }, { type: 'check', k: 'glu', label: 'Glu site' }];
    const fd = (v) => (v == null ? '' : (v >= 0 ? '+' : '−') + Math.abs(v / 1e6).toFixed(4));
    const len = (i) => ID.e[i] - ID.s[i] + 1;
    await tableView({ path: 'idrs', params, title: 'IDRs', unit: 'IDRs', items, n: ID.p.length, sortKey: 'ncpr', sortDir: 1, caption: 'IDRs', tsv: 'idr_explorer_idrs.tsv',
      lede: 'Every disordered region of ' + OV.settings.min_idr + ' or more residues (MobiDB), with how it is charged, how its charges are arranged, its sequence elements and its modification sites. Pick IDRs on the plot or with the filters; click one to open its protein there.' + docLink('charge'),
      filter: (i, st) => {
        if (st.minlen && len(i) < +st.minlen) return false;
        if (ps && !ps.has(ID.p[i])) return false;
        if (st.q && !qMatch(st.q)(ID.p[i])) return false;
        if (st.term && ID.term[i] !== st.term) return false;
        if (st.ct && ID.ct[i] !== +st.ct) return false;
        if (st.as && !ID.nas[i]) return false;
        if (st.rg && !ID.nrg[i]) return false;
        if (st.bp && !ID.nbp[i]) return false;
        if (st.mr && !ID.mR[i]) return false;
        if (st.glu && !ID.glu[i]) return false;
        if (st.ck2 && !ID.ck2[i]) return false;
        if (st.blocky && !(ID.dz[i] >= OV.settings.blocky_z * 100)) return false;
        if (st.fl && ID.zb[i] == null) return false;
        return true;
      },
      open: (i) => { location.hash = PR.acc[ID.p[i]] + '?f=' + ID.s[i] + '-' + ID.e[i]; },
      plot: {
        sub: 'Each IDR a point: grey all of them, coloured those your filters keep. Below the red line an IDR is acidic (net charge ≤ −0.10 per residue), above the blue one basic; near the corner few residues are charged. Colour by Das–Pappu class to see the diagram of states. Hover for the IDR, click to open it, drag to select IDRs for the table.',
        ax: 'fm', ay: 'fp', col: 'ct',
        axes: { fm: pctAxis('% acidic residues (D, E)', (i) => ID.fde[i] / 10), fp: pctAxis('% basic residues (K, R)', (i) => ID.frk[i] / 10), ncpr: { label: 'Net charge per residue', get: (i) => ID.ncpr[i] / 1000 },
          len: { label: 'Length (residues)', get: len, log: true }, kap: { label: 'Charge segregation (κ)', get: (i) => (ID.kap[i] == null ? null : ID.kap[i] / 1000) },
          dz: { label: 'Charge blocks (δ z)', get: (i) => z(ID.dz[i]) }, sites: { label: 'Modification sites in it', get: (i) => siteTotal(ID, i) },
          ph: { label: 'Phosphosites in it', get: (i) => ID.ph[i] }, mR: { label: 'Methylarginines in it', get: (i) => ID.mR[i] } },
        colours: { ct: { label: 'charge', cats: [['Acidic', '--acidic'], ['Basic', '--basic'], ['Mixed charge', '--hybrid'], ['Few charges', '--ink-2']], cat: (i) => [3, 3, 2, 0, 1][ID.ct[i]] },
          dp: { label: 'Das–Pappu class', cats: [1, 2, 3, 4, 5].map((r) => [REGION_SHORT[r], '--dp' + r]), cat: (i) => (ID.reg[i] || 1) - 1 },
          top: { label: 'acidic or basic flank', cats: [['At an acidic flank', '--acidic'], ['At a basic flank', '--basic'], ['At both', '--hybrid'], ['Neither', '--ink-2']], cat: (i) => (ID.top[i] === 3 ? 2 : ID.top[i] & 1 ? 0 : ID.top[i] & 2 ? 1 : 3) } },
        guides: chargeGuides,
        tip: (i) => '<b>' + esc(geneOf(ID.p[i])) + ' ' + ID.s[i] + '–' + ID.e[i] + '</b> (' + len(i) + ' aa)<br>' + CT[ID.ct[i]] + ': D/E ' + Math.round(ID.fde[i] / 10) + '%, K/R ' + Math.round(ID.frk[i] / 10) + '%, net ' + fsign(ID.ncpr[i] / 1000, 2) + ' per residue' +
          '<br>Das–Pappu ' + REGION_SHORT[ID.reg[i] || 1] + ' (' + REGION_SHAPE[ID.reg[i] || 1] + ')' +
          (ID.kap[i] != null ? '<br>charge segregation κ ' + f3(ID.kap[i]) : '') + (siteTotal(ID, i) ? '<br>' + fmt(siteTotal(ID, i)) + ' modification sites' : ''),
      },
      cols: [
        { k: 'gene', label: 'Protein', get: (i) => geneOf(ID.p[i]), td: (i) => geneCell(ID.p[i]) },
        { k: 'reg', label: 'IDR', tdCls: 'nw', get: (i) => ID.s[i], td: (i) => ID.s[i] + '–' + ID.e[i] + ' <span class="muted">' + ({ N: 'N-term', C: 'C-term', I: '', W: 'whole' }[ID.term[i]] || '') + '</span>' },
        { k: 'len', label: 'Length', num: true, get: len },
        { k: 'ncpr', label: 'Charge', num: true, get: (i) => ID.ncpr[i] / 1000, td: (i) => chargeCell(ID.fde[i] / 1000, ID.frk[i] / 1000, ID.ncpr[i] / 1000, ID.ct[i]), tip: 'D/E (red) and K/R (blue) share of residues, and net charge per residue; sorts by net charge' },
        { k: 'dp', label: 'Das–Pappu class', opt: true, on: true, get: (i) => ID.reg[i] || null, td: (i) => dpChip(ID.reg[i]), tip: 'Das & Pappu’s diagram-of-states class and the shape it predicts (hover a class)', menu: 'Das–Pappu class (predicted shape)' },
        { k: 'net', label: 'Net charge', num: true, get: (i) => Math.round(ID.ncpr[i] / 1000 * len(i)), td: (i) => fsign(Math.round(ID.ncpr[i] / 1000 * len(i)), 0), tip: 'K + R − D − E over the whole IDR' },
        { k: 'kap', label: 'Charge segregation (κ)', num: true, get: (i) => (ID.kap[i] == null ? null : ID.kap[i] / 1000), td: (i) => f3(ID.kap[i]), tip: 'κ (localCIDER): 0 positive and negative residues well mixed, 1 fully segregated; needs two charges of each sign' },
        { k: 'dz', label: 'Charge blocks (δ z)', num: true, opt: true, get: (i) => z(ID.dz[i]), td: (i) => fz(ID.dz[i]), tip: 'Charge patterning against 40 shuffles of the same residues; ≥ ' + OV.settings.blocky_z + ': blocks' },
        { k: 'om', label: 'Ω', num: true, opt: true, get: (i) => (ID.om[i] == null ? null : ID.om[i] / 1000), td: (i) => f3(ID.om[i]), tip: 'Ω (Martin et al. 2016): κ for charges and prolines together' },
        { k: 'hyd', label: 'Hydropathy', num: true, opt: true, get: (i) => (ID.hyd[i] == null ? null : ID.hyd[i] / 1000), td: (i) => f3(ID.hyd[i]), tip: 'Mean Kyte–Doolittle hydropathy, 0–1 (localCIDER)' },
        { k: 'dn', label: 'Δ net charge, modified', num: true, opt: true, tip: 'Net charge per residue with every observed site modified, minus unmodified (idrkit.ptm charges)', get: (i) => (ID.nmod[i] ? ID.dn[i] / 1000 : null), td: (i) => (ID.nmod[i] ? f3(ID.dn[i]) : '') },
        { k: 'dd', label: 'Δ δ, modified', num: true, opt: true, tip: 'δ with every observed site modified, minus unmodified: above 0, the modifications sharpen the charge blocks', get: (i) => (ID.nmod[i] ? ID.dd[i] / 1e6 : null), td: (i) => (ID.nmod[i] ? fd(ID.dd[i]) : '') },
        { k: 'lde', label: 'Longest D/E run', num: true, opt: true, get: (i) => ID.lde[i] },
        { k: 'nas', label: 'Acidic stretches', num: true, opt: true, get: (i) => ID.nas[i] }, { k: 'nbp', label: 'Basic patches', num: true, opt: true, get: (i) => ID.nbp[i] },
        { k: 'nrg', label: 'RG/RGG clusters', num: true, opt: true, get: (i) => ID.nrg[i] },
        { k: 'za', label: 'Acidic flank', num: true, opt: true, get: (i) => z(ID.za[i]), td: (i) => zChip(ID.za[i], 'acidic') + TOPB(ID.top[i] & 1), tip: 'Best acidic score among the domain flanks it overlaps' },
        { k: 'zb', label: 'Basic flank', num: true, opt: true, get: (i) => z(ID.zb[i]), td: (i) => zChip(ID.zb[i], 'basic') + TOPB(ID.top[i] & 2), tip: 'Best basic score among the domain flanks it overlaps' },
        { k: 'sites', label: 'Sites', num: true, get: (i) => siteTotal(ID, i), tip: 'Modification sites in the IDR (Columns: by kind)' },
      ].concat(KCOLS(ID, ['ph', 'mR'])),
    });
  }

  // Every modification site, as the Methylarginine Site Explorer's Sites table: where it sits and how it is charged there.
  async function sitesView(params) {
    view().innerHTML = '<p class="loading">Loading 334k sites…</p>';
    await Promise.all([loadProteins(), loadStats(), computeSelection()]);
    const ps = inProtSel(), keep = R.sel && R.sel.keep;
    const ctxBit = (c) => 1 << CTX_ORDER().indexOf(c);
    const CTX_SHOW = [['acidic_flank', 'acidic flank'], ['basic_flank', 'basic flank'], ['acidic_stretch', 'acidic stretch'], ['basic_patch', 'basic patch'], ['rg_cluster', 'RG/RGG cluster'], ['domain', 'domain']];
    const items = [{ type: 'text', k: 'q', label: 'Gene or accession', ph: 'e.g. HMGB1' },
      { type: 'select', k: 'kind', label: 'Modification', opts: OV.kinds.map(([k, l], b) => [b, l]) },
      { type: 'select', k: 'ct', label: 'In an IDR that is', opts: CT_OPTS.concat([[9, 'any IDR'], [0, 'not in an IDR']]) },
      { type: 'select', k: 'cx', label: 'Inside', opts: CTX_SHOW }];
    const kOf = (i) => OV.kinds[SI.k[i]];
    const KIN = OV.kinase_rows;
    await tableView({ path: 'sites', params, title: 'Sites', unit: 'modification sites', items, n: SI.p.length, sortKey: 'gene', sortDir: 1, caption: 'Sites', tsv: 'idr_explorer_sites.tsv',
      lede: 'Every modification site, with the charge of the disordered region it sits in and the sequence elements around it. Click a site to open it on its protein.' + docLink('kinds'),
      filter: (i, st) => {
        if (ps && !ps.has(SI.p[i])) return false;
        if (keep && !keep(i)) return false;
        if (st.q && !qMatch(st.q)(SI.p[i])) return false;
        if (st.kind !== '' && SI.k[i] !== +st.kind) return false;
        if (st.ct !== '') { const t = +st.ct; if (t === 9 ? !SI.ct[i] : SI.ct[i] !== t) return false; }
        if (st.cx && !(SI.ctx[i] & ctxBit(st.cx))) return false;
        return true;
      },
      open: (i) => { location.hash = PR.acc[SI.p[i]] + '?site=' + SI.pos[i]; },
      cols: [
        { k: 'gene', label: 'Protein', get: (i) => geneOf(SI.p[i]), td: (i) => geneCell(SI.p[i]) },
        { k: 'pos', label: 'Site', num: true, get: (i) => SI.pos[i], td: (i) => '<b>' + SI.r[i] + SI.pos[i] + '</b>' },
        { k: 'kind', label: 'Modification', get: (i) => kOf(i)[1], td: (i) => '<i class="kdot" style="background:var(--' + KCOL[kOf(i)[0]] + ')"></i>' + esc(kOf(i)[1]) },
        { k: 'ct', label: 'Its IDR', get: (i) => SI.ct[i], td: (i) => (SI.ct[i] ? ctChip(SI.ct[i]) : '<span class="muted">ordered</span>'), tip: 'The charge of the disordered region the site sits in' },
        { k: 'cx', label: 'Inside', get: (i) => CTX_SHOW.filter(([c]) => SI.ctx[i] & ctxBit(c)).map(([, l]) => l).join(', ') || null, td: (i) => CTX_SHOW.filter(([c]) => SI.ctx[i] & ctxBit(c)).map(([, l]) => '<span class="tag">' + esc(l) + '</span>').join(' ') },
        { k: 'km', label: 'Kinase motifs', get: (i) => (SI.km[i] > 0 ? KIN.filter((_, b) => SI.km[i] & (1 << b)).map((r) => r[1]).join(', ') : null), tip: 'Phosphosites: the representative kinases whose motif the site matches (Johnson et al. 2023, percentile ≥ ' + OV.kinase_match + ')' },
        { k: 'src', label: 'Reported by', get: (i) => srcNames(SI.src[i]) + (SI.conf[i] ? ', PTMeXchange ' + CONF[SI.conf[i]] : '') },
      ],
    });
  }

  // ---------------------------------------------------------------- protein page
  // What a reader of one protein wants to know, in order: which of its disordered regions are charged and how unusual
  // they are (highlights), where everything lies (the map, with a read-out of every residue under the pointer), each
  // charged segment's sequence and what its modifications do to its charge (segment cards, with a what-if), the sites
  // and their evidence, its structure, and its sequence. Built as the Methylarginine Site Explorer builds its page
  // (its header, facts, folds, marks, UniProt and AlphaFold loaders), with that page one click away.
  const KCOL = { ph: 'phospho', ck2: 'ck2', py: 'phospho', mR: 'meth', cit: 'cit', kac: 'kac', kacyl: 'kac', ub: 'ub', sumo: 'sumo', kme: 'kme', adpr: 'adpr', glu: 'glu' };
  // the light theme's colours (idr.css), for exports, which are always light
  const TRACK_LIGHT = { basic: '#2F6FD6', acidic: '#D64545', hybrid: '#7B3294', domain: '#B9DDE6', 'domain-ink': '#0F3B47', dis: '#8796A0', rg: '#F5D6DE', 'rg-ink': '#7A2E45', meth: '#B4485F', kme: '#A21CAF',
    phospho: '#6B7F88', ck2: '#B7791F', kac: '#EB6834', sumo: '#4D7C0F', glu: '#1B8A5A', cit: '#0E7490', adpr: '#78350F', ub: '#3F4C55', muted: '#5F6F78', 'ink-2': '#374151', rna: '#7C3AED', dna: '#B7791F', region: '#E7ECEF' };
  const tc = (name, exp) => (exp ? TRACK_LIGHT[name] || '#888' : C.cssVar('--' + name) || TRACK_LIGHT[name]);
  const CT_COL = ['', 'ink-2', 'hybrid', 'acidic', 'basic'];
  const MAP_KINDS = [['mR', 'Methyl-R'], ['cit', 'Citrulline'], ['kac', 'Acetyl-K'], ['kacyl', 'Acidic acyl-K'], ['ub', 'Ubiquitin-K'], ['sumo', 'SUMO-K'], ['kme', 'Methyl-K'],
    ['ph', 'Phospho-S/T'], ['py', 'Phospho-Y'], ['adpr', 'ADP-ribosyl'], ['glu', 'Glutamyl']];
  const MAP_DEFAULT = ['mR', 'cit', 'kac', 'kacyl', 'sumo', 'kme', 'ph', 'py', 'adpr', 'glu'];
  const mapKinds = () => { const s = store.get('map-kinds'); try { return new Set(s ? JSON.parse(s) : MAP_DEFAULT); } catch (e) { return new Set(MAP_DEFAULT); } };
  const MAP_ROWS = [['uni', 'Regions (UniProt)'], ['charge', 'Net charge'], ['idr', 'IDRs by charge'], ['flank', 'Acidic and basic flanks'], ['el', 'Sequence elements'], ['rna', 'RNA binding in cells'], ['dna', 'DNA binding in cells'], ['plddt', 'AlphaFold confidence'], ['ck2m', 'CK2 motifs without a site']];
  const mapRows = () => { const s = store.get('map-rows'); try { return new Set(s ? JSON.parse(s) : MAP_ROWS.map((r) => r[0])); } catch (e) { return new Set(MAP_ROWS.map((r) => r[0])); } };
  const PV = { p: null, zoom: null, labels: 'auto', pick: null, hl: null };
  const residueName = (seq, pos) => seq[pos - 1] + pos;
  // the charge a modification leaves on its residue (idrkit.ptm MODIFIED_CHARGE); kinds without one change no charge here
  const MOD_Q = { phospho: { S: -2, T: -2, Y: -2 }, acetyl: { K: 0 }, acyl_acidic: { K: -1 }, citrulline: { R: 0 }, adp_ribosyl: { S: -2, D: -2, E: -2 }, glutamyl: { E: -2 } };
  const RES_Q = { K: 1, R: 1, D: -1, E: -1 };
  const kindMod = (k) => (OV.kinds.find((x) => x[0] === k) || [])[3];
  const modQ = (k, res) => { const t = MOD_Q[kindMod(k)]; return t && res in t ? t[res] : null; };
  const SRC_LINK = { uniprot: (acc) => 'https://www.uniprot.org/uniprotkb/' + acc + '/entry#ptm_processing', iptmnet: (acc) => 'https://research.bioinformatics.udel.edu/iptmnet/entry/' + acc + '/',
    release: (acc) => PTM_SITE + '#' + acc, ptmx: () => 'https://www.ebi.ac.uk/pride/archive/projects/PXD079218', rebak2024: () => 'https://doi.org/10.1038/s41594-024-01214-9',
    hendriks2021: () => 'https://doi.org/10.1038/s41467-021-26172-4', koenig2026: () => 'https://doi.org/10.1038/s41594-026-01877-6', ttll4: () => 'https://www.ebi.ac.uk/pride/archive/projects/PXD065251' };
  const srcLinks = (bits, conf, acc) => OV.sources.filter((s) => bits & s[1]).map(([k, , l]) => '<a href="' + SRC_LINK[k](acc) + '" target="_blank" rel="noopener">' + esc(l) + '</a>' + (k === 'ptmx' && conf ? ' ' + CONF[conf] : '')).join(', ');
  const srcText = (v) => srcNames(v[1]) + (v[2] ? ', PTMeXchange ' + CONF[v[2]] : '');
  const PLDDT = [[90, '#0053D6', 'very high (≥ 90)'], [70, '#65CBF3', 'confident (70–90)'], [50, '#FFDB13', 'low (50–70)'], [0, '#FF7D45', 'very low (< 50)']];
  const plddtCol = (v) => PLDDT.find((b) => v >= b[0])[1];
  const XL_PEP = ['RBDmap fragment', 'CAPRI peptide', 'pCLAP peptide', 'RBR-ID+ stretch'];
  // ranks against every IDR: the share of IDRs more acidic (basic) than this one
  let RANK = null;
  const rankOf = (i) => {
    if (!RANK) { const v = Int32Array.from(ID.ncpr).sort(); RANK = v; }
    const x = ID.ncpr[i]; let lo = 0, hi = RANK.length;
    while (lo < hi) { const m = (lo + hi) >> 1; if (RANK[m] < x) lo = m + 1; else hi = m; }
    let hi2 = lo; while (hi2 < RANK.length && RANK[hi2] === x) hi2++;
    return { below: lo / RANK.length, above: (RANK.length - hi2) / RANK.length };   // shares more acidic, more basic
  };
  const rankText = (i) => {
    const r = rankOf(i), t = ID.ct[i];
    if (t === 3) return 'more acidic than ' + (100 * (1 - r.below)).toFixed(r.below < 0.01 ? 1 : 0) + '% of IDRs';
    if (t === 4) return 'more basic than ' + (100 * (1 - r.above)).toFixed(r.above < 0.01 ? 1 : 0) + '% of IDRs';
    return '';
  };
  const netOf = (i) => Math.round(ID.ncpr[i] / 1000 * (ID.e[i] - ID.s[i] + 1));

  // the charged segments of a protein: its IDRs, and its acidic and basic flanks (the top tenth by either score) where
  // they are not already most of an IDR, as a short acidic tail beside a fold is not an IDR of 20 residues (MYB 1-34);
  // each with its charge from its sequence and what touches it
  function segmentsOf(ctx) {
    const { Rp, S, idrs, names, seq } = ctx, xl = ctx.xl || {};
    const dn = (d) => (names[d[2]] ? names[d[2]][1] : 'a domain');
    const make = (o) => {
      const { s, e } = o, sub = seq.slice(s - 1, e), len = e - s + 1;
      const fm = sub.replace(/[^DE]/g, '').length / len, fp = sub.replace(/[^KR]/g, '').length / len, net = Math.round((fp - fm) * len);
      const sites = [];
      for (const [k] of MAP_KINDS) for (const v of S[k] || []) if (v[0] >= s && v[0] <= e) sites.push({ k, v, pos: v[0] });
      sites.sort((x, y) => x.pos - y.pos);
      const before = (Rp.dom || []).filter((d) => d[0] > e && d[0] - e <= 10).map((d) => 'before ' + dn(d) + ' (' + d[0] + '–' + d[1] + ')');
      const after = (Rp.dom || []).filter((d) => d[1] < s && s - d[1] <= 10).map((d) => 'after ' + dn(d) + ' (' + d[0] + '–' + d[1] + ')');
      return Object.assign({ len, fm, fp, net, ncpr: net / len, ct: ctOf(fm + fp, net / len), sites, next: after.concat(before),
        rna: (xl.rna || []).filter((q) => q >= s && q <= e).length, dna: (xl.dna || []).filter((q) => q >= s && q <= e).length }, o);
    };
    const out = idrs.map((i) => make({ id: 'i' + i, i, s: ID.s[i], e: ID.e[i], kind: 'IDR', fl: (Rp.fl || []).filter((f) => Math.min(ID.e[i], f[2]) - Math.max(ID.s[i], f[1]) + 1 >= 10) }));
    // acidic and basic flanks: of overlapping windows of the same kind, the best-scoring one (a union would dilute it)
    for (const [bit, what] of [[1, 'Acidic flank'], [2, 'Basic flank']]) {
      const sc = (f) => (bit === 1 ? f[9] : f[8]) || 0;
      const spans = (Rp.fl || []).filter((f) => (f[10] || 0) & bit).sort((x, y) => x[1] - y[1]), groups = [];
      for (const f of spans) { const last = groups[groups.length - 1]; if (last && f[1] <= last.end) { last.fs.push(f); last.end = Math.max(last.end, f[2]); } else groups.push({ end: f[2], fs: [f] }); }
      const joined = groups.map((gr) => { const best = gr.fs.slice().sort((x, y) => sc(y) - sc(x))[0]; return [best[1], best[2], [best]]; });
      for (const [a, b, fs] of joined) {
        const inIdr = idrs.some((i) => Math.min(b, ID.e[i]) - Math.max(a, ID.s[i]) + 1 >= 0.5 * (b - a + 1));
        if (!inIdr) out.push(make({ id: 'f' + a + '-' + b, i: null, s: a, e: b, kind: what, bit, fl: fs }));
      }
    }
    return out.sort((x, y) => x.s - y.s);
  }

  // ---- the map
  function protMapSvg(ctx, o = {}) {
    const { p, seq, Rp, S, idrs, names } = ctx, n = seq.length, exp = !!o.exp, rowsOn = mapRows();
    const pal = exp ? C.EXPORT_PAL : Object.assign(C.screenPalette(), { bg: C.cssVar('--surface') });
    const z = PV.zoom || [1, n], v0 = Math.max(1, z[0]), v1 = Math.min(n, Math.max(z[1], v0 + 10));
    const W = o.W || 1000, ML = Math.min(150, Math.round(W * 0.18)), MR = 16;
    const X = (pos) => ML + ((pos - v0) / (v1 - v0 + 1)) * (W - ML - MR);
    const cx = (pos) => Math.max(ML, Math.min(W - MR, X(pos)));
    const vis = (a, b) => b >= v0 && a <= v1;
    const span = (a, b) => (a === b ? '' + a : a + '–' + b);
    const tip = (html) => (exp ? '' : ' data-tip="' + esc(html) + '"');
    const legend = [];
    const keyOnce = (k, item) => { if (!legend.some((x) => x[2] === k)) legend.push(item.concat(k)); };
    let g = '';
    // sites: one head per kind at each residue, in levels where they would touch
    const on = mapKinds(), ck2 = new Set((S.ck2 || []).map((v) => v[0]));
    const sites = [];
    for (const [k] of MAP_KINDS) if (on.has(k)) for (const v of S[k] || []) if (v[0] >= v0 && v[0] <= v1) sites.push({ k: k === 'ph' && ck2.has(v[0]) ? 'ck2' : k, kind: k, v, pos: v[0] });
    sites.sort((a, b) => a.pos - b.pos || MAP_KINDS.findIndex((x) => x[0] === a.kind) - MAP_KINDS.findIndex((x) => x[0] === b.kind));
    const R0 = 4.6, levels = [];
    const heads = sites.map((s) => {
      const hx = X(s.pos + 0.5);
      let lv = 0;
      while (levels[lv] != null && hx - levels[lv] < 2 * R0 + 1 && lv < 11) lv++;
      levels[lv] = hx;
      return Object.assign(s, { hx, lv });
    });
    const nLv = Math.max(1, levels.length);
    const LFS = 9.5, LINES = PV.labels === 'all' ? 10 : 3, placed = [];
    if (PV.labels !== 'none' && heads.length && (PV.labels === 'all' || heads.length <= 80)) {
      const ends = new Array(LINES).fill(-1e9), byPos = new Map();
      for (const h of heads) { if (!byPos.has(h.pos)) byPos.set(h.pos, []); byPos.get(h.pos).push(h); }
      for (const [pos, hs] of byPos) {
        const text = residueName(seq, pos), w = C.textW(text, LFS, pal), lx = Math.max(ML + w / 2, Math.min(W - MR - w / 2, hs[0].hx));
        const line = ends.findIndex((e) => lx - w / 2 > e + 4);
        if (line < 0) continue;
        ends[line] = lx + w / 2;
        placed.push({ hs, text, lx, line });
      }
    }
    const nLines = placed.length ? Math.max(...placed.map((q) => q.line)) + 1 : 0;
    const y0 = exp ? 64 : 6, headTop = y0 + (nLines ? nLines * 12 + 8 : 0), base = headTop + 14 + nLv * 11 + 12;
    if (exp) {
      const title = geneOf(p) + ' · ' + PR.name[p], sub = PR.acc[p] + ' · ' + fmt(n) + ' aa' + (PV.zoom ? ' · residues ' + v0 + '–' + v1 : '') + ' · IDR Explorer';
      g += '<text x="0" y="24" font-size="' + C.fitFs(title, 18, W - MR, pal, 700) + '" font-weight="700" fill="' + pal.ink + '">' + esc(title) + '</text><text x="0" y="44" font-size="11.5" fill="' + pal.muted + '">' + esc(sub) + '</text>';
    }
    const rowLab = (txt, yy, help) => { g += '<text x="0" y="' + (yy + 4) + '" font-size="' + C.fitFs(txt, 10, ML - 8, pal, 600) + '" font-weight="600" fill="' + pal.muted + '"' + (help && !exp ? ' class="row-lab" data-tip="' + esc('<b>' + esc(txt) + '</b><br>' + help) + '"' : '') + '>' + esc(txt) + '</text>'; };
    if (heads.length) rowLab('Sites', headTop + 8, 'Modification sites, one head per kind; levels stack heads that would touch. Click a head for the site.');
    const dis = { list: Rp.dis || [], source: 'MobiDB' };
    const chain = C.chainSvg(pal, X, v0, v1, base, 8, dis, (svg, a, b) => (exp ? svg : svg.replace(/\/>$/, ' class="hit" data-seg="' + a + '-' + b + '" data-what="Disordered"/>')));
    g += chain.bars + chain.waves;
    rowLab('Chain', base, 'Bar: ordered; wave: disordered (MobiDB); boxes: InterPro domains and repeats.');
    if (dis.list.some(([a, b]) => vis(a, b))) keyOnce('dis', [C.SWATCH.disorder(pal), 'Disordered (MobiDB)']);
    const below = [];
    for (const d of (Rp.dom || []).filter((d) => vis(d[0], d[1]))) {
      const a = cx(d[0]), w = Math.max(2, cx(d[1] + 1) - a), nm = names[d[2]] ? names[d[2]][1] : '';
      g += C.MARK.domain(pal, a, base - 9, w, 18, false, exp ? '' : ' class="hit" data-seg="' + d[0] + '-' + d[1] + '" data-what="' + esc(nm) + '"');
      const lab = nm.replace(/ domain$/i, '');
      if (lab && C.textW(lab, 10, pal, 600) < w - 6) g += '<text x="' + (a + w / 2).toFixed(1) + '" y="' + (base + 3.5) + '" font-size="10" font-weight="600" text-anchor="middle" fill="' + pal.domainInk + '" pointer-events="none">' + esc(lab) + '</text>';
      else if (lab) below.push([a + w / 2, lab]);
    }
    if ((Rp.dom || []).some((d) => vis(d[0], d[1]))) keyOnce('dom', [C.SWATCH.domain(pal), 'Domain or repeat (InterPro)']);
    const headY = (h) => headTop + 14 + (nLv - 1 - h.lv) * 11;
    for (const h of heads) {
      const hy = headY(h), col = tc(KCOL[h.k], exp), sel = PV.pick && PV.pick.pos === h.pos && !exp;
      g += '<line x1="' + h.hx.toFixed(1) + '" x2="' + h.hx.toFixed(1) + '" y1="' + (base - 4) + '" y2="' + hy + '" stroke="' + pal.stem + '" stroke-width="1"/>';
      g += C.MARK.head(pal, h.hx, hy, R0, col, false, exp ? '' : ' class="site-head" data-pos="' + h.pos + '"', sel ? pal.ink : null);
    }
    for (const q of placed) {
      const ly = y0 + 9 + q.line * 12, h = q.hs[0];
      g += '<line x1="' + h.hx.toFixed(1) + '" y1="' + (headY(h) - R0 - 1) + '" x2="' + q.lx.toFixed(1) + '" y2="' + (ly + 3) + '" stroke="' + pal.stem + '" stroke-width=".7"/>';
      g += '<text x="' + q.lx.toFixed(1) + '" y="' + ly + '" font-size="' + LFS + '" text-anchor="middle" fill="' + pal.ink + '">' + esc(q.text) + '</text>';
    }
    const kindsDrawn = new Set(heads.map((h) => h.k));
    for (const k of ['mR', 'cit', 'kac', 'ub', 'sumo', 'kme', 'ph', 'ck2', 'py', 'adpr', 'glu']) if (kindsDrawn.has(k) || (k === 'kac' && kindsDrawn.has('kacyl'))) {
      keyOnce('k-' + k, [C.SWATCH.head(pal, tc(KCOL[k], exp), false), { mR: 'Methyl-R', cit: 'Citrulline', kac: 'Acetyl or acyl-K', ub: 'Ubiquitin-K', sumo: 'SUMO-K', kme: 'Methyl-K', ph: 'Phospho-S/T', ck2: 'Phospho-S/T, CK2-attributed', py: 'Phospho-Y', adpr: 'ADP-ribosyl', glu: 'Glutamyl' }[k]]);
    }
    let y = base + 22;
    if (below.length) {
      let right = -1e9;
      for (const [mx, text] of below.sort((a, b) => a[0] - b[0])) { const w = C.textW(text, 9.5, pal); if (mx - w / 2 > right + 6) { g += '<text x="' + mx.toFixed(1) + '" y="' + (y - 4) + '" font-size="9.5" text-anchor="middle" fill="' + pal.muted + '">' + esc(text) + '</text>'; right = mx + w / 2; } }
      y += 12;
    }
    const grid = (yy) => { g += '<line x1="' + ML + '" x2="' + (W - MR) + '" y1="' + yy + '" y2="' + yy + '" stroke="' + pal.grid + '"/>'; };
    const box = (a, b, yy, h, fill, attr, extra) => '<rect x="' + cx(a).toFixed(1) + '" y="' + (yy - h / 2).toFixed(1) + '" width="' + Math.max(1.5, cx(b + 1) - cx(a)).toFixed(1) + '" height="' + h + '" rx="2" fill="' + fill + '"' + (extra || '') + (exp ? '' : attr) + '/>';
    const tick = (pos, yy, h, fill, attr) => '<rect x="' + (X(pos + 0.5) - 1).toFixed(1) + '" y="' + (yy - h / 2) + '" width="2" height="' + h + '" fill="' + fill + '"' + (exp ? '' : attr || '') + '/>';
    const label = (a, b, yy, text, fill) => { const w = cx(b + 1) - cx(a); if (text && C.textW(text, 9, pal) < w - 6) g += '<text x="' + (cx(a) + 3).toFixed(1) + '" y="' + (yy + 3.2) + '" font-size="9" fill="' + fill + '" pointer-events="none">' + esc(text) + '</text>'; };
    // UniProt regions, motifs and compositional bias (read live)
    const uni = ctx.uni && ctx.uni.feats ? ctx.uni.feats.filter((f) => ['Region', 'Motif', 'Compositional bias'].includes(f.type) && f.desc !== 'Disordered') : [];
    if (rowsOn.has('uni') && uni.length) {
      // packed onto up to four lanes; a label sits inside its box, else after it, else (no room) only in the tooltip
      const items = uni.filter((f) => vis(f.start, f.end)).sort((a, b) => a.start - b.start || b.end - a.end), lanes = [];
      for (const f of items) {
        const a = cx(f.start), b = cx(f.end + 1), tw = C.textW(f.desc, 9, pal) + 8, inside = tw < b - a;
        const reach = (lab) => (lab && !inside ? b + tw : b);
        let k = lanes.findIndex((r) => r < a - 2);
        if (k < 0 && lanes.length < 4) { lanes.push(-1e9); k = lanes.length - 1; }
        f.lab = k >= 0 && (inside || b + tw < W - MR);
        if (k < 0) { k = lanes.indexOf(Math.min(...lanes)); f.lab = false; }
        lanes[k] = Math.max(lanes[k], reach(f.lab));
        f.lane = k;
      }
      rowLab('Regions', y, 'Regions, motifs and compositional bias from the UniProt entry, read live.');
      for (const f of items) {
        const yy = y + f.lane * 13;
        g += box(f.start, f.end, yy, 10, f.type === 'Motif' ? pal.bg : tc('region', exp), ' class="hit" data-seg="' + f.start + '-' + f.end + '" data-what="' + esc(f.desc) + '"', f.type === 'Motif' ? ' stroke="' + pal.muted + '"' : '');
        const w = cx(f.end + 1) - cx(f.start);
        if (f.lab) g += '<text x="' + (cx(f.start) + (C.textW(f.desc, 9, pal) + 8 < w ? 3 : w + 3)).toFixed(1) + '" y="' + (yy + 3.2) + '" font-size="9" fill="' + pal.ink + '" pointer-events="none">' + esc(f.desc) + '</text>';
      }
      if (items.length) keyOnce('uni', [C.SWATCH.fill(tc('region', exp)), 'Region (UniProt)']);
      y += 20 + Math.max(0, lanes.length - 1) * 13;
    }
    // net charge: 5-residue mean, up basic, down acidic
    if (rowsOn.has('charge')) {
      const Hc = 34, mid = y + 4, bins = Math.min(v1 - v0 + 1, 500), bw = (W - ML - MR) / bins;
      rowLab('Net charge', mid, 'Net charge of the 5 residues around each position (K, R +1; D, E −1; H +0.5): up basic, down acidic.');
      g += '<line x1="' + ML + '" x2="' + (W - MR) + '" y1="' + mid + '" y2="' + mid + '" stroke="' + pal.grid + '"/>';
      chargeProfile(seq.slice(v0 - 1, v1), bins).forEach((v, b) => { if (!v) return; const hh = Math.abs(v) * Hc / 2; g += '<rect x="' + (ML + b * bw).toFixed(1) + '" y="' + (v > 0 ? mid - hh : mid).toFixed(1) + '" width="' + (bw + 0.2).toFixed(2) + '" height="' + hh.toFixed(1) + '" fill="' + tc(v > 0 ? 'basic' : 'acidic', exp) + '"/>'; });
      y += Hc / 2 + 18;
    }
    // IDRs coloured by charge, labelled with their net charge
    if (rowsOn.has('idr') && idrs.length) {
      rowLab('IDRs', y, 'Each disordered region (≥ ' + OV.settings.min_idr + ' aa) coloured by its charge: acidic (net charge ≤ −0.10 per residue), basic (≥ +0.10), mixed (a quarter or more charged, balanced), or few charges; labelled with its net charge. Click one for its sequence.');
      grid(y);
      for (const i of idrs.filter((i) => vis(ID.s[i], ID.e[i]))) {
        const t = ID.ct[i];
        g += box(ID.s[i], ID.e[i], y, 13, tc(CT_COL[t], exp), ' class="hit" data-seg="' + ID.s[i] + '-' + ID.e[i] + '" data-idr="' + i + '" data-what="IDR"', t === 1 ? ' fill-opacity="0.35"' : '');
        label(ID.s[i], ID.e[i], y, fsign(netOf(i), 0), t === 1 ? pal.ink : '#fff');
      }
      for (const t of [3, 4, 2, 1]) if (idrs.some((i) => ID.ct[i] === t)) keyOnce('ct' + t, [C.SWATCH.fill(tc(CT_COL[t], exp), t === 1 ? 0.35 : null), 'IDR, ' + CT[t].toLowerCase()]);
      y += 20;
    }
    // acidic and basic flanks: the top tenth by each score, solid
    const tops = (Rp.fl || []).filter((f) => (f[10] || 0) && vis(f[1], f[2]));
    if (rowsOn.has('flank') && (Rp.fl || []).some((f) => f[10])) {
      rowLab('Flanks', y, 'Domain flanks in the top tenth of all flanks by the acidic (red) or basic (blue) score: the acidic and basic flanks. Click one for its sequence.');
      grid(y);
      for (const f of tops) for (const [bit, col, what] of [[1, 'acidic', 'Acidic flank'], [2, 'basic', 'Basic flank']]) if (f[10] & bit) {
        g += box(f[1], f[2], y + (f[10] === 3 ? (bit === 1 ? -3 : 3) : 0), f[10] === 3 ? 6 : 10, tc(col, exp), ' class="hit" data-seg="' + f[1] + '-' + f[2] + '" data-what="' + what + ' (' + (f[0] ? 'after' : 'before') + ' the domain at ' + f[3] + '–' + f[4] + ')"');
        keyOnce('fl' + bit, [C.SWATCH.fill(tc(col, exp)), what]);
      }
      y += 18;
    }
    // sequence elements on one row
    const els = [['as', 'acidic', 'Acidic stretch'], ['bp', 'basic', 'Basic patch'], ['rg', 'rg', 'RG/RGG cluster']];
    if (rowsOn.has('el') && els.some(([k]) => (Rp[k] || []).length)) {
      rowLab('Elements', y, 'Acidic stretches (red), basic patches (blue: net charge +3 or more in 10 residues) and RG/RGG clusters (pink), from idrkit.');
      grid(y);
      for (const [k, col, what] of els) for (const d of (Rp[k] || []).filter((d) => vis(d[0], d[1]))) {
        g += box(d[0], d[1], y, k === 'rg' ? 6 : 10, tc(col, exp), ' class="hit" data-seg="' + d[0] + '-' + d[1] + '" data-what="' + what + '"', k === 'rg' ? ' stroke="' + tc('rg-ink', exp) + '" stroke-width="0.6"' : ' fill-opacity="0.6"');
        keyOnce('el-' + k, [C.SWATCH.fill(tc(col, exp), k === 'rg' ? null : 0.6), what]);
      }
      y += 18;
    }
    // RNA and DNA binding in living cells (the Methylarginine Site Explorer's compilation)
    const xl = ctx.xl || {};
    if (rowsOn.has('rna') && ((xl.rna || []).length || (xl.rpep || []).length)) {
      rowLab('RNA in cells', y, 'Residues UV-crosslinked to RNA in living cells (ticks: Kramer 2014, XRNAX, RBS-ID, pRBS-ID, pepR-MS compendium) and RNA-bound stretches (boxes: RBDmap fragments, CAPRI and pCLAP peptides, RBR-ID+ stretches), as the Methylarginine Site Explorer compiles them.');
      grid(y);
      for (const [a, b, k] of (xl.rpep || []).filter((x) => vis(x[0], x[1]))) g += box(a, b, y, 10, tc('rna', exp), ' class="hit" data-seg="' + a + '-' + b + '" data-what="RNA-bound ' + XL_PEP[k] + '"', ' fill-opacity="0.22"');
      for (const q of (xl.rna || []).filter((q) => q >= v0 && q <= v1)) g += tick(q, y, 14, tc('rna', exp));
      keyOnce('rna', [C.SWATCH.tick(tc('rna', exp)), 'RNA crosslink in cells (box: RNA-bound stretch)']);
      y += 18;
    }
    if (rowsOn.has('dna') && ((xl.dna || []).length || (xl.dpep || []).length || (xl.db || []).length)) {
      rowLab('DNA in cells', y, 'Residues and peptides crosslinked to DNA in cells (Trendel et al. 2025), and UniProt DNA-binding regions (outlined).');
      grid(y);
      for (const [a, b] of (xl.db || []).filter((x) => vis(x[0], x[1]))) g += box(a, b, y, 12, 'none', ' class="hit" data-seg="' + a + '-' + b + '" data-what="DNA-binding region (UniProt)"', ' stroke="' + tc('dna', exp) + '" stroke-width="1.2"');
      for (const [a, b] of (xl.dpep || []).filter((x) => vis(x[0], x[1]))) g += box(a, b, y, 10, tc('dna', exp), ' class="hit" data-seg="' + a + '-' + b + '" data-what="DNA-crosslinked peptide"', ' fill-opacity="0.3"');
      for (const q of (xl.dna || []).filter((q) => q >= v0 && q <= v1)) g += tick(q, y, 14, tc('dna', exp));
      keyOnce('dna', [C.SWATCH.tick(tc('dna', exp)), 'DNA crosslink in cells (outline: UniProt DNA-binding region)']);
      y += 18;
    }
    // AlphaFold confidence along the chain (read live)
    if (rowsOn.has('plddt') && ctx.af && ctx.af.plddt && ctx.af.plddt.size) {
      rowLab('AlphaFold', y, 'AlphaFold DB model confidence (pLDDT) along the chain, read live: blue confident, yellow and orange low (often disorder).');
      const bins = Math.min(v1 - v0 + 1, 600), bw = (W - ML - MR) / bins;
      for (let b = 0; b < bins; b++) { const pos = Math.round(v0 + (b + 0.5) * (v1 - v0 + 1) / bins), v = ctx.af.plddt.get(pos); if (v != null) g += '<rect x="' + (ML + b * bw).toFixed(1) + '" y="' + (y - 5) + '" width="' + (bw + 0.3).toFixed(2) + '" height="10" fill="' + plddtCol(v) + '"/>'; }
      keyOnce('af', [(x2, y2) => PLDDT.map((b, k) => '<rect x="' + (x2 + k * 5) + '" y="' + (y2 - 5) + '" width="5" height="10" fill="' + b[1] + '"/>').join(''), 'AlphaFold pLDDT (blue high, orange low)']);
      y += 18;
    }
    if (rowsOn.has('ck2m') && (S.ck2m || []).length) {
      rowLab('CK2 motifs', y, 'Serines and threonines CK2 scores at the 99th percentile or above (Johnson et al. 2023) with no phosphosite on record: candidates.');
      grid(y);
      for (const v of S.ck2m.filter((v) => v[0] >= v0 && v[0] <= v1)) g += C.MARK.head(pal, X(v[0] + 0.5), y, 3.6, tc('ck2', exp), true, exp ? '' : ' class="site-head" data-pos="' + v[0] + '"');
      keyOnce('ck2m', [C.SWATCH.head(pal, tc('ck2', exp), true), 'Strong CK2 motif, no site reported']);
      y += 18;
    }
    const hl = PV.hl;
    if (hl && vis(hl[0], hl[1]) && !exp) g = '<rect x="' + cx(hl[0]).toFixed(1) + '" y="' + headTop + '" width="' + Math.max(2, cx(hl[1] + 1) - cx(hl[0])).toFixed(1) + '" height="' + (y - headTop) + '" fill="' + pal.accent + '" fill-opacity="0.08" stroke="' + pal.accent + '" stroke-opacity="0.45"/>' + g;
    const step = [5, 10, 20, 50, 100, 200, 250, 500, 1000, 2000, 5000].find((s) => (v1 - v0) / s <= 10) || 10000;
    g += '<line x1="' + ML + '" x2="' + (W - MR) + '" y1="' + y + '" y2="' + y + '" stroke="' + pal.muted + '"/>';
    for (let t = Math.ceil(v0 / step) * step; t <= v1; t += step) g += '<line x1="' + X(t + 0.5).toFixed(1) + '" x2="' + X(t + 0.5).toFixed(1) + '" y1="' + y + '" y2="' + (y + 4) + '" stroke="' + pal.muted + '"/><text x="' + X(t + 0.5).toFixed(1) + '" y="' + (y + 15) + '" font-size="10" text-anchor="middle" fill="' + pal.muted + '">' + t + '</text>';
    y += 24;
    let H = y;
    if (exp) { const lg = C.svgLegend(pal, legend.map(([d, l]) => [d, l]), 0, y + 4, W - MR, { fs: 10 }); g += lg.svg; H = Math.ceil(y + lg.h + 14); }
    else g += '<line id="xhair" x1="0" x2="0" y1="' + headTop + '" y2="' + (y - 24) + '" stroke="' + pal.ink + '" stroke-opacity="0.35" visibility="hidden" pointer-events="none"/>';
    return { W, H, v0, v1, x0: ML, x1: W - MR, len: n, legend,
      svg: '<svg xmlns="http://www.w3.org/2000/svg" class="pmap" width="' + W + '" height="' + H + '" viewBox="0 0 ' + W + ' ' + H + '" font-family=\'' + (pal.font || C.EXPORT_PAL.font).replace(/'/g, '"') + '\' role="img" aria-label="' + esc('Protein map of ' + geneOf(p)) + '">' +
        (exp ? '<rect width="' + W + '" height="' + H + '" fill="' + pal.bg + '"/>' : '') + g + '</svg>' };
  }
  const EMBOSS = { K: 1, R: 1, H: 0.5, D: -1, E: -1 };
  function chargeProfile(seq, bins) {
    const n = seq.length, w = 2, q = Array.from(seq, (c) => EMBOSS[c] || 0), m = new Float32Array(n);
    for (let i = 0; i < n; i++) { let s = 0, k = 0; for (let j = Math.max(0, i - w); j <= Math.min(n - 1, i + w); j++) { s += q[j]; k++; } m[i] = s / k; }
    const out = [];
    for (let b = 0; b < bins; b++) { const a = Math.floor(b * n / bins), e = Math.max(a + 1, Math.floor((b + 1) * n / bins)); let s = 0; for (let i = a; i < e; i++) s += m[i]; out.push(s / (e - a)); }
    return out;
  }
  // everything recorded at one residue: the read-out under the pointer
  function residueReadout(ctx, pos) {
    const { seq, Rp, S, idrs, names } = ctx, out = ['<b>' + residueName(seq, pos) + '</b>'];
    const i = idrs.find((k) => ID.s[k] <= pos && ID.e[k] >= pos);
    if (i != null) out.push('IDR ' + ID.s[i] + '–' + ID.e[i] + ': ' + CT[ID.ct[i]].toLowerCase() + ', net ' + fsign(netOf(i), 0) + ', Das–Pappu ' + REGION_SHORT[ID.reg[i] || 1]);
    else if ((Rp.dis || []).some(([a, b]) => a <= pos && b >= pos)) out.push('Disordered (MobiDB), in a stretch under ' + OV.settings.min_idr + ' residues');
    for (const d of Rp.dom || []) if (d[0] <= pos && d[1] >= pos) out.push(esc(names[d[2]] ? names[d[2]][1] : 'Domain') + ' ' + d[0] + '–' + d[1]);
    for (const f of Rp.fl || []) if ((f[10] || 0) && f[1] <= pos && f[2] >= pos) out.push(((f[10] & 1) ? 'Acidic' : '') + (f[10] === 3 ? ' and basic' : (f[10] & 2) ? 'Basic' : '') + ' flank ' + f[1] + '–' + f[2]);
    for (const [k, l] of [['as', 'Acidic stretch'], ['bp', 'Basic patch'], ['rg', 'RG/RGG cluster']]) for (const d of Rp[k] || []) if (d[0] <= pos && d[1] >= pos) out.push(l + ' ' + d[0] + '–' + d[1]);
    if (ctx.uni && ctx.uni.feats) { const u = ctx.uni.feats.filter((f) => ['Region', 'Motif', 'Compositional bias'].includes(f.type) && f.desc !== 'Disordered' && f.start <= pos && f.end >= pos).map((f) => f.desc); if (u.length) out.push('UniProt: ' + esc(u.slice(0, 3).join(' · '))); }
    const here = [];
    for (const [k, l] of MAP_KINDS) for (const v of S[k] || []) if (v[0] === pos) here.push(esc(l) + ' <span class="muted">(' + esc(srcText(v)) + ')</span>');
    if (here.length) out.push(here.join('<br>'));
    const xl = ctx.xl || {};
    if ((xl.rna || []).includes(pos)) out.push('RNA crosslink in cells');
    if ((xl.dna || []).includes(pos)) out.push('DNA crosslink in cells');
    if (ctx.af && ctx.af.plddt && ctx.af.plddt.get(pos) != null) out.push('AlphaFold pLDDT ' + ctx.af.plddt.get(pos).toFixed(0));
    let ch = 0;
    for (let k = Math.max(1, pos - 2); k <= Math.min(seq.length, pos + 2); k++) ch += RES_Q[seq[k - 1]] || 0;
    out.push('<span class="muted">Net charge of residues ' + Math.max(1, pos - 2) + '–' + Math.min(seq.length, pos + 2) + ': ' + fsign(ch, 0) + '</span>');
    return out.join('<br>');
  }
  // a stretch coloured by charge, modified residues underlined; mods: Map(pos -> new charge) shows a what-if
  function seqHtml(ctx, s, e, o = {}) {
    const marks = new Map();
    for (const k of Object.keys(ctx.S)) if (k !== 'ck2m' && k !== 'ck2') for (const v of ctx.S[k]) { if (!marks.has(v[0])) marks.set(v[0], []); marks.get(v[0]).push(k); }
    let html = '';
    for (let i = s; i <= e; i++) {
      const c = ctx.seq[i - 1], q = o.mods && o.mods.has(i) ? o.mods.get(i) : RES_Q[c] || 0;
      const cls = (q > 0 ? 'b' : q < 0 ? 'a' : '') + (marks.has(i) ? ' s' : '') + (o.mods && o.mods.has(i) ? ' m' : '') + (i === o.pick ? ' pick' : '');
      if ((i - s) % 60 === 0) html += (i > s ? '<br>' : '') + '<span class="n">' + String(i).padStart(5, ' ') + '</span> ';
      html += cls.trim() || marks.has(i) ? '<span class="' + cls.trim() + '" data-tip="' + esc(c + i + (marks.has(i) ? ': ' + marks.get(i).map(kindLabel).join(', ') : '') + (o.mods && o.mods.has(i) ? ' · charge ' + fsign(o.mods.get(i), 0) + ' when modified' : '')) + '">' + c + '</span>' : c;
      if ((i - s + 1) % 10 === 0 && (i - s + 1) % 60 !== 0 && i < e) html += ' ';
    }
    return '<div class="seq">' + html + '</div>';
  }
  // a net-charge profile of a segment, before (grey) and after (colour) its checked modifications
  function segProfile(ctx, sg, mods) {
    const W = 600, H = 70, s = sg.s, e = sg.e, n = e - s + 1;
    const q0 = [], q1 = [];
    for (let i = s; i <= e; i++) { const c = ctx.seq[i - 1], a = RES_Q[c] || 0; q0.push(a); q1.push(mods.has(i) ? mods.get(i) : a); }
    const sm = (q) => q.map((_, i) => { let t = 0, k = 0; for (let j = Math.max(0, i - 2); j <= Math.min(n - 1, i + 2); j++) { t += q[j]; k++; } return t / k; });
    const a = sm(q0), b = sm(q1), X = (i) => 4 + (i + 0.5) / n * (W - 8), Y = (v) => H / 2 - v * (H / 2 - 4) / 2;
    const path = (v) => v.map((x, i) => (i ? 'L' : 'M') + X(i).toFixed(1) + ' ' + Y(x).toFixed(1)).join('');
    const area = (v, sign) => 'M' + X(0).toFixed(1) + ' ' + Y(0) + v.map((x, i) => 'L' + X(i).toFixed(1) + ' ' + Y(sign * x > 0 ? x : 0).toFixed(1)).join('') + 'L' + X(n - 1).toFixed(1) + ' ' + Y(0) + 'Z';
    return '<svg class="segprof" viewBox="0 0 ' + W + ' ' + H + '" preserveAspectRatio="none" role="img" aria-label="Net charge along the segment">' +
      '<line x1="0" x2="' + W + '" y1="' + Y(0) + '" y2="' + Y(0) + '" class="sp-zero"/><path d="' + area(b, 1) + '" class="sp-b"/><path d="' + area(b, -1) + '" class="sp-a"/>' +
      (mods.size ? '<path d="' + path(a) + '" class="sp-before"/>' : '') + '</svg>';
  }
  function segCard(ctx, sg) {
    const i = sg.i, isIdr = i != null;
    const changeable = sg.sites.filter((x) => modQ(x.k, ctx.seq[x.pos - 1]) != null);
    const fl = isIdr ? sg.fl.filter((f) => f[10]).map((f) => ((f[10] & 1) ? 'acidic' : '') + (f[10] === 3 ? ' and basic' : (f[10] & 2) ? 'basic' : '') + ' flank ' + f[1] + '–' + f[2]) : [];
    const score = !isIdr ? sg.kind.split(' ')[0].toLowerCase() + ' score ' + pctTop(Math.max(...sg.fl.map((f) => (sg.bit === 1 ? f[9] : f[8]) || 0))) + ' of all flanks' : '';
    const kinds = {};
    for (const x of sg.sites) kinds[x.k] = (kinds[x.k] || 0) + 1;
    return '<div class="seg-card' + (isIdr ? '' : ' seg-flank') + '" id="seg-' + sg.id + '"><div class="seg-head"><div><h4>' + ctChip(sg.ct) + ' ' + (isIdr ? dpChip(ID.reg[i]) + ' ' : '') + esc(sg.kind) + ' ' + sg.s + '–' + sg.e + ' <span class="muted">· ' + sg.len + ' aa</span></h4>' +
      '<p class="sub">' + chargeCell(sg.fm, sg.fp, sg.ncpr, sg.ct) + ' net charge <b>' + fsign(sg.net, 0) + '</b>' + (isIdr && rankText(i) ? ' · ' + rankText(i) : '') + (score ? ' · ' + score : '') +
      (isIdr && ID.kap[i] != null ? ' · charge segregation κ ' + f3(ID.kap[i]) : '') + (sg.next.length ? ' · ' + esc(sg.next.join('; ')) : '') + (fl.length ? ' · ' + esc(fl.join('; ')) : '') +
      (!isIdr ? ' · <span class="muted">a domain flank, not an IDR of ' + OV.settings.min_idr + ' or more residues</span>' : '') +
      (sg.rna ? ' · ' + sg.rna + ' RNA-crosslinked residue' + (sg.rna > 1 ? 's' : '') : '') + (sg.dna ? ' · ' + sg.dna + ' DNA-crosslinked' : '') + '</p></div>' +
      '<div class="exp"><button type="button" data-segzoom="' + sg.s + '-' + sg.e + '">Show on the map</button><button type="button" data-segcopy="' + sg.s + '-' + sg.e + '">Copy sequence</button></div></div>' +
      '<div class="seg-seq" data-seq>' + seqHtml(ctx, sg.s, sg.e) + '</div>' +
      '<div class="seg-prof" data-prof>' + segProfile(ctx, sg, new Map()) + '</div>' +
      (sg.sites.length ? '<div class="seg-mods"><span class="flabel">Its ' + sg.sites.length + ' modification site' + (sg.sites.length > 1 ? 's' : '') + '</span> ' +
        Object.entries(kinds).map(([k, c]) => '<i class="kdot" style="background:var(--' + KCOL[k] + ')"></i>' + c + ' ' + esc(kindLabel(k))).join(' · ') +
        (changeable.length ? '<div class="whatif"><span class="flabel">What if they were modified?</span> <button class="btn quiet small" type="button" data-modall>All ' + changeable.length + '</button><button class="btn quiet small" type="button" data-modnone>None</button><span class="wi-out" data-out></span><div class="wi-list">' +
          changeable.map((x) => '<label class="check" data-tip="' + esc(kindLabel(x.k) + ' at ' + residueName(ctx.seq, x.pos) + ': charge ' + fsign(RES_Q[ctx.seq[x.pos - 1]] || 0, 0) + ' → ' + fsign(modQ(x.k, ctx.seq[x.pos - 1]), 0) + ' (' + srcText(x.v) + ')') + '"><input type="checkbox" data-mod="' + x.pos + ':' + x.k + '"> <i class="kdot" style="background:var(--' + KCOL[x.k] + ')"></i>' + residueName(ctx.seq, x.pos) + '</label>').join('') + '</div></div>' :
          '<p class="note">None of them changes charge in this model (methylation keeps it; ubiquitin and SUMO attach a protein).</p>') + '</div>' : '<p class="note">No modification site reported here.</p>') + '</div>';
  }
  function bindSegCard(ctx, el, sg) {
    const out = $('[data-out]', el);
    const update = () => {
      const mods = new Map();
      for (const cb of $$('[data-mod]', el)) if (cb.checked) { const [pos, k] = cb.dataset.mod.split(':'); mods.set(+pos, modQ(k, ctx.seq[+pos - 1])); }
      let net1 = 0, charged = 0;
      for (let k = sg.s; k <= sg.e; k++) { const q = mods.has(k) ? mods.get(k) : RES_Q[ctx.seq[k - 1]] || 0; net1 += q; if (q) charged++; }
      const ncpr1 = net1 / sg.len, t1 = ctOf(charged / sg.len, ncpr1);
      if (out) out.innerHTML = mods.size ? 'net charge <b>' + fsign(sg.net, 0) + ' → ' + fsign(net1, 0) + '</b> (' + fsign(sg.ncpr, 2) + ' → ' + fsign(ncpr1, 2) + ' per residue)' + (t1 !== sg.ct ? ': ' + CT[sg.ct].toLowerCase() + ' → <b>' + CT[t1].toLowerCase() + '</b>' : ', still ' + CT[t1].toLowerCase()) : '';
      $('[data-seq]', el).innerHTML = seqHtml(ctx, sg.s, sg.e, { mods });
      $('[data-prof]', el).innerHTML = segProfile(ctx, sg, mods);
    };
    $$('[data-mod]', el).forEach((cb) => cb.addEventListener('change', update));
    const all = $('[data-modall]', el), none = $('[data-modnone]', el);
    if (all) all.addEventListener('click', () => { $$('[data-mod]', el).forEach((cb) => { cb.checked = true; }); update(); });
    if (none) none.addEventListener('click', () => { $$('[data-mod]', el).forEach((cb) => { cb.checked = false; }); update(); });
    $('[data-segcopy]', el).addEventListener('click', (e) => { try { navigator.clipboard.writeText(ctx.seq.slice(sg.s - 1, sg.e)); e.currentTarget.textContent = 'Copied'; } catch (err) { /* no clipboard */ } });
  }
  // the panel under the map: a site (its kinds, sources, motifs, occupancy and where it sits) or a stretch
  function pickPanel(ctx) {
    const box = $('#site-card'), pk = PV.pick;
    if (!box) return;
    if (!pk) { box.hidden = true; return; }
    box.hidden = false;
    const { seq, S, Rp, idrs } = ctx, n = seq.length;
    const close = '<button type="button" id="pick-close">Close</button>';
    const zoomBtn = (a, b) => '<button type="button" data-zoom="' + Math.max(1, a - 15) + '-' + Math.min(n, b + 15) + '">Zoom the map to it</button>';
    if (pk.pos != null) {
      const pos = pk.pos, rows = [];
      for (const [k, l] of MAP_KINDS.concat([['ck2', 'CK2-attributed phospho']])) for (const v of S[k] || []) if (v[0] === pos) rows.push([k, l, v]);
      const dark = (S.ck2m || []).find((v) => v[0] === pos), kin = rows.find(([k, , v]) => k === 'ph' && v.length > 4), KIN = OV.kinase_rows;
      const i = idrs.find((k) => ID.s[k] <= pos && ID.e[k] >= pos), res = seq[pos - 1];
      const dq = rows.map(([k]) => modQ(k, res)).find((q) => q != null);
      box.innerHTML = '<div class="card-head"><div><h3>' + residueName(seq, pos) + '</h3><p class="sub">' + (rows.length ? rows.filter(([k]) => k !== 'ck2').map(([, l]) => esc(l)).join(' · ') : 'no site reported') + (dark ? ' · strong CK2 motif' : '') + '</p></div><div class="exp">' + zoomBtn(pos, pos) + close + '</div></div>' +
        '<dl class="kv">' + rows.filter(([k]) => k !== 'ck2').map(([k, l, v]) => '<dt>' + esc(l) + '</dt><dd>' + srcLinks(v[1], v[2], PR.acc[ctx.p]) + (k === 'kac' && typeof v[3] === 'number' ? '; stoichiometry ' + (100 * v[3]).toFixed(3) + '% (Hansen 2019)' : '') + (typeof v[3] === 'string' && v[3] ? '; ' + esc(v[3]) : '') + '</dd>').join('') +
        (dq != null ? '<dt>Charge</dt><dd>' + res + ' ' + fsign(RES_Q[res] || 0, 0) + ' → ' + fsign(dq, 0) + ' when modified' + (i != null ? '; its IDR’s net charge ' + fsign(netOf(i), 0) + ' → ' + fsign(netOf(i) + dq - (RES_Q[res] || 0), 0) : '') + '</dd>' : '') +
        (rows.some(([k]) => k === 'ck2') ? '<dt>Kinase</dt><dd>CK2, as UniProt reports it</dd>' : '') +
        (kin ? '<dt>Kinase motifs</dt><dd>CK2 percentile ' + kin[2][3] + '; matches (≥ ' + OV.kinase_match + '): ' + (KIN.filter((_, b) => kin[2][4] & (1 << b)).map((r) => r[1]).join(', ') || 'none') + '; best-scoring kinase ' + esc(B.top_kinases[kin[2][5]]) + ' (' + kin[2][6] + ')</dd>' : '') +
        (kin && kin[2][7] != null ? '<dt>Occupancy</dt><dd>' + (kin[2][7] / 10).toFixed(1) + '% (median of 33 cell lines, Koenig 2026)</dd>' : '') +
        (dark ? '<dt>CK2 motif</dt><dd>percentile ' + dark[1] + ', no phosphosite on record</dd>' : '') +
        '<dt>Where</dt><dd>' + residueReadout(ctx, pos).split('<br>').slice(1, -1).join('<br>') + '</dd></dl>' +
        seqHtml(ctx, Math.max(1, pos - 15), Math.min(n, pos + 15), { pick: pos }) +
        (rows.some(([k]) => k === 'mR') ? '<p class="note"><a href="' + PTM_SITE + '#' + PR.acc[ctx.p] + '" target="_blank" rel="noopener">Its methylarginine evidence in the Methylarginine Site Explorer ↗</a></p>' : '');
    } else {
      const [a, b] = pk.seg, sub = seq.slice(a - 1, b), q = Array.from(sub, (c) => RES_Q[c] || 0), net = q.reduce((s, v) => s + v, 0);
      const fm = sub.replace(/[^DE]/g, '').length / sub.length, fp = sub.replace(/[^KR]/g, '').length / sub.length;
      box.innerHTML = '<div class="card-head"><div><h3>' + esc(pk.what || 'Residues') + ' ' + a + '–' + b + '</h3><p class="sub">' + (b - a + 1) + ' residues · ' + chargeCell(fm, fp, net / sub.length, ctOf(fm + fp, net / sub.length)) + ' net charge <b>' + fsign(net, 0) + '</b>' +
        (pk.idr != null ? ' · <a href="#" data-goseg="' + pk.idr + '">its segment card</a>' : '') + '</p></div><div class="exp">' + zoomBtn(a, b) + close + '</div></div>' + seqHtml(ctx, a, b);
      const gs = $('[data-goseg]', box);
      if (gs) gs.addEventListener('click', (e) => { e.preventDefault(); openSeg('i' + gs.dataset.goseg); });
    }
    $('#pick-close', box).addEventListener('click', () => { PV.pick = null; drawProtMap(ctx); pickPanel(ctx); });
    $$('[data-zoom]', box).forEach((bt) => bt.addEventListener('click', () => { const [x, y] = bt.dataset.zoom.split('-').map(Number); PV.zoom = [x, y]; drawProtMap(ctx); $('#pmap-card').scrollIntoView({ block: 'nearest' }); }));
  }
  const openSeg = (id) => {
    const fold = $('details[data-fold="p-seg"]');
    if (fold) fold.open = true;
    const more = $('#seg-more');
    if (!$('#seg-' + id) && more) more.click();
    const el = $('#seg-' + id);
    if (el) { el.scrollIntoView({ block: 'start', behavior: 'smooth' }); el.classList.remove('flash'); void el.offsetWidth; el.classList.add('flash'); }
  };
  function drawProtMap(ctx) {
    const box = $('#map-box');
    if (!box) return;
    const w = box.clientWidth, out = protMapSvg(ctx, { W: w > 0 && w < 760 ? Math.max(560, Math.round(w)) : 1000 });
    box.innerHTML = out.svg;
    const pal = C.screenPalette();
    $('#fig-legend').innerHTML = out.legend.map(([draw, label]) => '<span>' + C.swatchHtml(draw, pal) + esc(label) + '</span>').join('');
    $('#zoom-reset').hidden = !PV.zoom;
    $('#map-sub').textContent = PR.acc[ctx.p] + ' · ' + fmt(out.len) + ' aa' + (PV.zoom ? ' · showing ' + out.v0 + '–' + out.v1 : '') + (ctx.uni ? ' · regions from UniProt' : '') + (ctx.af ? ' · confidence from AlphaFold DB' : '');
    const svg = $('svg', box), xh = $('#xhair', svg);
    const toX = (e) => { const r = svg.getBoundingClientRect(); return ((e.clientX - r.left) / r.width) * out.W; };
    const toPos = (vx) => Math.floor(out.v0 + ((vx - out.x0) / (out.x1 - out.x0)) * (out.v1 - out.v0 + 1));
    let start = null, brush = null, moved = false;
    svg.addEventListener('pointerdown', (e) => {
      if (e.button !== 0) return;
      const vx = toX(e);
      if (vx < out.x0 || vx > out.x1) return;
      start = vx; moved = false;
      dragging(true);
      svg.setPointerCapture(e.pointerId);
      brush = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
      brush.setAttribute('y', 0); brush.setAttribute('height', out.H); brush.setAttribute('x', vx); brush.setAttribute('width', 0);
      brush.setAttribute('fill', C.cssVar('--accent')); brush.setAttribute('fill-opacity', '0.12'); brush.setAttribute('pointer-events', 'none');
      svg.appendChild(brush);
    });
    svg.addEventListener('pointermove', (e) => {
      const vx = toX(e);
      if (start != null) {
        const cur = Math.max(out.x0, Math.min(out.x1, vx));
        if (Math.abs(cur - start) > 4) { moved = true; C.hideTip(); }
        brush.setAttribute('x', Math.min(start, cur)); brush.setAttribute('width', Math.abs(cur - start));
        return;
      }
      if (vx < out.x0 || vx > out.x1) { xh.setAttribute('visibility', 'hidden'); C.hideTip(); return; }
      const pos = Math.max(out.v0, Math.min(out.v1, toPos(vx))), px = out.x0 + ((pos - out.v0 + 0.5) / (out.v1 - out.v0 + 1)) * (out.x1 - out.x0);
      xh.setAttribute('x1', px); xh.setAttribute('x2', px); xh.setAttribute('visibility', 'visible');
      const t = e.target.closest('.site-head');
      C.showTip(e, residueReadout(ctx, t ? +t.dataset.pos : pos) + '<br><span class="muted">' + (t ? 'Click for the site' : e.target.closest('[data-seg]') ? 'Click to read this stretch' : 'Drag to zoom') + '</span>');
    });
    svg.addEventListener('pointerleave', () => { xh.setAttribute('visibility', 'hidden'); C.hideTip(); });
    svg.addEventListener('pointerup', (e) => {
      dragging(false);
      if (start == null) return;
      const cur = Math.max(out.x0, Math.min(out.x1, toX(e)));
      const a = toPos(Math.min(start, cur)), b = toPos(Math.max(start, cur));
      start = null;
      brush.remove();
      if (moved && b - a >= 8) { PV.zoom = [a, b]; drawProtMap(ctx); return; }
      const el = document.elementFromPoint(e.clientX, e.clientY), head = el && el.closest('.site-head'), seg = el && el.closest('[data-seg]');
      if (head) { PV.pick = { pos: +head.dataset.pos }; C.hideTip(); drawProtMap(ctx); pickPanel(ctx); $('#site-card').scrollIntoView({ block: 'nearest' }); }
      else if (seg) { const [x, y] = seg.dataset.seg.split('-').map(Number); PV.pick = { seg: [x, y], what: seg.dataset.what, idr: seg.dataset.idr != null ? +seg.dataset.idr : null }; C.hideTip(); pickPanel(ctx); $('#site-card').scrollIntoView({ block: 'nearest' }); }
    });
    svg.addEventListener('dblclick', () => { if (PV.zoom) { PV.zoom = null; drawProtMap(ctx); } });
  }

  // ---- the 3D structure: AlphaFold DB model (3Dmol.js, as the Methylarginine Site Explorer loads it)
  function structureView(ctx, mode) {
    const host = $('#af-viewer'), note = $('#af-note');
    if (!host) return;
    if (!ctx.af) { host.innerHTML = '<p class="note">' + (ctx.afTried ? 'No AlphaFold DB model for ' + esc(PR.acc[ctx.p]) + '.' : 'Looking up AlphaFold DB…') + '</p>'; return; }
    C.load3Dmol().then((lib) => C.fetchText(ctx.af.entry.pdbUrl).then((pdb) => {
      if (!$('#af-viewer')) return;
      host.innerHTML = '';
      const v = lib.createViewer(host, { backgroundColor: C.cssVar('--surface') || 'white' });
      v.addModel(pdb, 'pdb');
      const colourOf = (atom) => {
        const pos = atom.resi, c = ctx.seq[pos - 1];
        if (mode === 'plddt') return plddtCol(atom.b);
        if (mode === 'idr') { const i = ctx.idrs.find((k) => ID.s[k] <= pos && ID.e[k] >= pos); return i == null ? '#C7D4D9' : TRACK_LIGHT[CT_COL[ID.ct[i]]]; }
        return c === 'D' || c === 'E' ? TRACK_LIGHT.acidic : c === 'K' || c === 'R' ? TRACK_LIGHT.basic : '#D9DEE2';
      };
      v.setStyle({}, { cartoon: { colorfunc: colourOf, thickness: 0.4 } });
      const shown = mapKinds(), ck2 = new Set((ctx.S.ck2 || []).map((x) => x[0]));
      for (const [k] of MAP_KINDS) if (shown.has(k)) for (const x of ctx.S[k] || []) v.addStyle({ resi: x[0], atom: 'CA' }, { sphere: { radius: 1.2, color: TRACK_LIGHT[KCOL[k === 'ph' && ck2.has(x[0]) ? 'ck2' : k]] } });
      v.setHoverable({}, true, (atom, viewer) => { if (!atom.label) atom.label = viewer.addLabel(ctx.seq[atom.resi - 1] + atom.resi + (atom.b != null ? ' · pLDDT ' + atom.b.toFixed(0) : ''), { position: atom, backgroundColor: '#111827', fontColor: 'white', fontSize: 12 }); },
        (atom, viewer) => { if (atom.label) { viewer.removeLabel(atom.label); delete atom.label; } });
      v.zoomTo(); v.zoom(1.4); v.render();
      note.textContent = 'AlphaFold DB ' + ctx.af.entry.entryId + (ctx.af.entry.latestVersion ? ' v' + ctx.af.entry.latestVersion : '') + ' · spheres: the modification sites shown on the map · hover a residue for its name and confidence.';
    })).catch(() => { host.innerHTML = '<p class="note">The structure viewer could not load.</p>'; });
  }

  async function proteinView(key, params, token) {
    await loadProteins();
    const K0 = key.toUpperCase(), p = byAcc.has(K0) ? byAcc.get(K0) : byGene.get(K0);
    if (p == null) { view().innerHTML = '<div class="panel">No reviewed human protein called <b>' + esc(key) + '</b>. <a href="#proteins?q=' + encodeURIComponent(key) + '">Search the protein table</a></div>'; return; }
    const acc = PR.acc[p];
    const [sh] = await Promise.all([loadShard(acc), loadIdrs()]);
    if (token !== routeToken) return;
    const Rp = sh.prot[acc], seq = Rp.seq, n = seq.length, S = Rp.s || {};
    const idrs = [];
    for (let i = 0; i < ID.p.length; i++) if (ID.p[i] === p) idrs.push(i);
    const ctx = { p, seq, Rp, S, idrs, names: sh.names, xl: Rp.xl || null, uni: null, af: null, afTried: false };
    if (PV.p !== p) Object.assign(PV, { p, zoom: null, pick: null });
    const hl = (params.get('f') || '').match(/^(\d+)-(\d+)$/), site = +(params.get('site') || 0);
    PV.hl = hl ? [+hl[1], +hl[2]] : null;
    if (PV.hl && n > 400) PV.zoom = [Math.max(1, PV.hl[0] - 60), Math.min(n, PV.hl[1] + 60)];
    if (PV.hl) PV.pick = { seg: PV.hl, what: 'Selected', idr: idrs.find((i) => ID.s[i] === PV.hl[0] && ID.e[i] === PV.hl[1]) ?? null };
    if (site) { PV.pick = { pos: site }; if (n > 300) PV.zoom = [Math.max(1, site - 60), Math.min(n, site + 60)]; }
    UI.detail = p;
    if (!UI.pOpen || UI.pOpen.p !== p) UI.pOpen = { p, open: new Map() };
    const segs = segmentsOf(ctx);
    // highlights: what is unusual about this protein, each opening its place
    const acidic = idrs.filter((i) => ID.ct[i] === 3).sort((a, b) => ID.ncpr[a] - ID.ncpr[b])[0], basic = idrs.filter((i) => ID.ct[i] === 4).sort((a, b) => ID.ncpr[b] - ID.ncpr[a])[0];
    const glyph = (i) => '<span class="hg-seq">' + Array.from(seq.slice(ID.s[i] - 1, ID.e[i]).match(/[DEKR]{2}.{0,4}/) ? seq.slice(ID.s[i] - 1, ID.e[i]).match(/[DEKR]{2}.{0,4}/)[0] : seq.slice(ID.s[i] - 1, ID.s[i] + 5)).slice(0, 6).map((c) => '<i class="' + ('DE'.includes(c) ? 'a' : 'KR'.includes(c) ? 'b' : '') + '">' + c + '</i>').join('') + '</span>';
    const hi = [];
    if (acidic != null) hi.push({ go: 'seg:i' + acidic, glyph: glyph(acidic), big: fsign(netOf(acidic), 0), dir: 'down', label: 'Most acidic IDR, ' + ID.s[acidic] + '–' + ID.e[acidic] + ' (' + (ID.e[acidic] - ID.s[acidic] + 1) + ' aa)', kind: rankText(acidic), tip: 'Net charge ' + fsign(netOf(acidic), 0) + ': D/E ' + Math.round(ID.fde[acidic] / 10) + '%, K/R ' + Math.round(ID.frk[acidic] / 10) + '%. Click for its sequence and modifications.' });
    if (basic != null) hi.push({ go: 'seg:i' + basic, glyph: glyph(basic), big: fsign(netOf(basic), 0), dir: 'up', label: 'Most basic IDR, ' + ID.s[basic] + '–' + ID.e[basic] + ' (' + (ID.e[basic] - ID.s[basic] + 1) + ' aa)', kind: rankText(basic), tip: 'Net charge ' + fsign(netOf(basic), 0) + ': K/R ' + Math.round(ID.frk[basic] / 10) + '%, D/E ' + Math.round(ID.fde[basic] / 10) + '%. Click for its sequence and modifications.' });
    const inSegs = new Map();
    for (const sg of segs) for (const x of sg.sites) inSegs.set(x.pos + x.k, x);
    const inIdr = [...inSegs.values()], dq = inIdr.reduce((t, x) => { const q = modQ(x.k, seq[x.pos - 1]); return q == null ? t : t + q - (RES_Q[seq[x.pos - 1]] || 0); }, 0);
    const nSites = OV.kinds.filter(([k]) => k !== 'ck2').reduce((s, [k]) => s + (S[k] || []).length, 0);
    if (nSites) hi.push({ go: 'sites', glyph: C.aGlyph.pins(), big: fmt(inIdr.length) + ' of ' + fmt(nSites), dir: 'up', label: 'modification sites in its IDRs and charged flanks', kind: dq ? 'all modified, they would change the charge by ' + fsign(dq, 0) : 'none changes charge in the model', tip: 'Sites inside disordered regions; the charge change counts phosphate −2, acetyl-K −1, acidic acyl-K −2, citrulline −1, ADP-ribosyl −2, glutamyl −1 relative to the unmodified residue.' });
    const topFl = (Rp.fl || []).filter((f) => f[10]).sort((a, b) => Math.max(b[8] || 0, b[9] || 0) - Math.max(a[8] || 0, a[9] || 0))[0];
    if (topFl) { const ac = (topFl[10] & 1) && (topFl[9] || 0) >= (topFl[8] || 0); hi.push({ go: 'fl:' + topFl[1] + '-' + topFl[2], glyph: C.aGlyph.struct(), big: pctTop(ac ? topFl[9] : topFl[8]), dir: ac ? 'down' : 'up', label: (ac ? 'Acidic' : 'Basic') + ' flank ' + topFl[1] + '–' + topFl[2], kind: (topFl[0] ? 'after' : 'before') + ' the domain at ' + topFl[3] + '–' + topFl[4], tip: 'Its best domain flank by the ' + (ac ? 'acidic' : 'basic') + ' score, against all 121,444 flanks.' }); }
    const xlr = ctx.xl ? (ctx.xl.rna || []).filter((q) => segs.some((sg) => q >= sg.s && q <= sg.e)).length : 0;
    if (ctx.xl && (ctx.xl.rna || []).length) hi.push({ go: 'map', glyph: C.aGlyph.cls(), big: fmt(xlr) + ' of ' + fmt(ctx.xl.rna.length), dir: 'up', label: 'RNA-crosslinked residues in its IDRs', kind: 'in living cells', tip: 'Residues UV-crosslinked to RNA in cells (the Methylarginine Site Explorer’s compilation); the map shows where.' });
    const ext = (href, label) => '<a class="btn quiet small" href="' + href + '" target="_blank" rel="noopener">' + label + ' ↗</a>';
    const prow = (k, title, inner, open, sum) => C.foldHtml(k, title, inner, open, sum);
    view().innerHTML = '<nav class="crumbs"><a href="#proteins">Proteins</a> › <b>' + esc(geneOf(p)) + '</b></nav>' +
      '<div class="detail-head"><div class="detail-title"><h2>' + esc(geneOf(p)) + '</h2><span class="pn">' + esc(PR.name[p]) + '</span></div>' +
      '<p class="pfunc" id="pfunc"><span class="muted">Function: loading from UniProt…</span></p>' +
      '<div class="detail-links"><span class="tag mono">' + esc(acc) + '</span><span class="tag">' + fmt(n) + ' aa</span><span class="tag">' + PR.dis[p] + '% disordered</span>' +
      OV.classes.filter((c, b) => PR.cls[p] & (1 << b)).map((c) => '<span class="tag">' + esc(c[1]) + '</span>').join('') +
      ext('https://www.uniprot.org/uniprotkb/' + acc + '/entry', 'UniProt') + '<a class="btn small" href="' + PTM_SITE + '#' + acc + '" target="_blank" rel="noopener">Methylarginine Site Explorer ↗</a>' +
      '<details class="links-more"><summary class="btn quiet small">More links</summary><div class="links-menu">' + ext('https://alphafold.ebi.ac.uk/entry/' + acc, 'AlphaFold') + ext('https://mobidb.org/' + acc, 'MobiDB') +
      ext('https://www.ebi.ac.uk/interpro/protein/UniProt/' + acc + '/', 'InterPro') + ext('https://research.bioinformatics.udel.edu/iptmnet/entry/' + acc + '/', 'iPTMnet') + '</div></details></div></div>' +
      (hi.length ? C.aHighlights(hi) : '') +
      prow('p-map', 'Protein map', '<div class="card fig-card" id="pmap-card"><div class="card-head"><div><h3>' + esc(geneOf(p)) + ': charge, disorder and modifications</h3><p class="sub" id="map-sub"></p></div>' +
        '<div class="exp" role="group" aria-label="Figure export"><button type="button" id="zoom-reset" hidden>Reset zoom</button><button type="button" id="fig-svg">SVG</button><button type="button" id="fig-png">PNG</button></div></div>' +
        '<div class="fig-controls"><div class="seg" role="group" aria-label="Site labels"><span>Site labels</span>' + [['auto', 'Auto'], ['all', 'All'], ['none', 'None']].map(([k, l]) => '<button type="button" data-lab="' + k + '" aria-pressed="' + (PV.labels === k) + '">' + l + '</button>').join('') + '</div>' +
        '<details class="col-menu"><summary class="btn quiet small">Sites</summary><div class="col-pop" role="group" aria-label="Kinds of site drawn">' +
        MAP_KINDS.map(([k, l]) => '<label class="check"><input type="checkbox" data-kind="' + k + '"' + (mapKinds().has(k) ? ' checked' : '') + ((S[k] || []).length ? '' : ' disabled') + '> ' + esc(l) + ' <span class="muted">' + fmt((S[k] || []).length) + '</span></label>').join('') + '</div></details>' +
        '<details class="col-menu"><summary class="btn quiet small">Rows</summary><div class="col-pop" role="group" aria-label="Rows drawn">' + MAP_ROWS.map(([k, l]) => '<label class="check"><input type="checkbox" data-row="' + k + '"' + (mapRows().has(k) ? ' checked' : '') + '> ' + esc(l) + '</label>').join('') + '</div></details></div>' +
        '<div class="map-box" id="map-box"></div><details class="key-fold" open><summary>Key</summary><div class="legend" id="fig-legend"></div></details>' +
        '<p class="note">Point anywhere for what is at that residue; click a site or a stretch to read it below; drag across to zoom, double-click to zoom out. SVG and PNG reproduce this view with a title and key.</p></div>', true, fmt(nSites) + ' sites · ' + idrs.length + ' IDRs') +
      '<div class="card" id="site-card" hidden></div>' +
      '<div class="p-folds">' +
      prow('p-seg', 'Charged segments', '<div class="card"><p class="sub">Each disordered region, and each acidic or basic flank outside one: its sequence coloured by charge (acidic red, basic blue, modified residues underlined), its net charge along the sequence (5-residue mean), and what its modifications would do to its charge. Tick sites to see.</p><div id="segs"></div></div>', true,
        [3, 4, 2].filter((t) => segs.some((x) => x.ct === t)).map((t) => segs.filter((x) => x.ct === t).length + ' ' + CT[t].toLowerCase()).join(' · ') || (segs.length ? 'few charges' : 'no IDR of ' + OV.settings.min_idr + ' aa or more')) +
      prow('p-sites', 'Modification sites', '<div id="psites"></div>', false, fmt(nSites) + ' sites') +
      prow('p-fl', 'Domain flanks', '<div id="pfl"></div>', false, fmt((Rp.fl || []).length) + ' flanks') +
      prow('p-3d', '3D structure', '<div class="card"><div class="card-head"><div><h3>AlphaFold model</h3><p class="sub" id="af-note">Loads when opened.</p></div><div class="seg" role="group" aria-label="Colour the structure by"><span>Colour by</span>' +
        [['charge', 'Charge'], ['idr', 'IDRs'], ['plddt', 'Confidence']].map(([k, l], j) => '<button type="button" data-afc="' + k + '" aria-pressed="' + (j === 0) + '">' + l + '</button>').join('') + '</div></div><div class="af-viewer" id="af-viewer"></div></div>', false, 'AlphaFold DB') +
      prow('p-seq', 'Sequence', '<div class="card"><div class="card-head"><div><h3>Sequence</h3><p class="sub">Basic residues blue, acidic red, modified residues underlined; hover a residue for its marks.</p></div><div class="exp"><button type="button" id="copy-fasta">Copy FASTA</button></div></div>' + seqHtml(ctx, 1, n) + '</div>', false, fmt(n) + ' aa') + '</div>';
    document.title = geneOf(p) + ' · IDR Explorer';
    drawProtMap(ctx);
    pickPanel(ctx);
    // segment cards: the most charged first when there are many
    const segOrder = segs.slice().sort((a, b) => a.s - b.s), MANY = 8;
    const strong = segs.length > MANY ? new Set(segs.slice().sort((a, b) => Math.abs(b.net) + b.sites.length - Math.abs(a.net) - a.sites.length).slice(0, MANY).map((x) => x.id)) : null;
    const drawSegs = (all) => {
      const list = all || !strong ? segOrder : segOrder.filter((x) => strong.has(x.id));
      $('#segs').innerHTML = list.map((sg) => segCard(ctx, sg)).join('') + (!all && strong ? '<button class="btn" type="button" id="seg-more">Show all ' + segs.length + ' segments (showing the ' + MANY + ' most charged or modified)</button>' : '') +
        (segs.length ? '' : '<p class="note">No disordered region of ' + OV.settings.min_idr + ' or more residues, and no acidic or basic flank.</p>');
      for (const sg of list) bindSegCard(ctx, $('#seg-' + sg.id), sg);
      $$('[data-segzoom]').forEach((b) => b.addEventListener('click', () => { const [a, c] = b.dataset.segzoom.split('-').map(Number); PV.zoom = [Math.max(1, a - 20), Math.min(n, c + 20)]; PV.hl = [a, c]; drawProtMap(ctx); $('#pmap-card').scrollIntoView({ block: 'start', behavior: 'smooth' }); }));
      const more = $('#seg-more');
      if (more) more.addEventListener('click', () => drawSegs(true));
    };
    drawSegs(false);
    $$('.a-hi-c[data-go]').forEach((b) => b.addEventListener('click', () => {
      const [k, v] = b.dataset.go.split(':');
      if (k === 'seg') { const i = +v.slice(1); PV.zoom = n > 300 ? [Math.max(1, ID.s[i] - 30), Math.min(n, ID.e[i] + 30)] : null; PV.hl = [ID.s[i], ID.e[i]]; drawProtMap(ctx); openSeg(v); }
      else if (k === 'fl') { const [a, c] = v.split('-').map(Number); PV.zoom = n > 300 ? [Math.max(1, a - 40), Math.min(n, c + 40)] : null; PV.hl = [a, c]; drawProtMap(ctx); const sg = segs.find((x) => x.i == null && x.s <= a && x.e >= c); if (sg) openSeg(sg.id); else $('#pmap-card').scrollIntoView({ block: 'start', behavior: 'smooth' }); }
      else if (k === 'sites') { const f = $('details[data-fold="p-sites"]'); f.open = true; f.scrollIntoView({ block: 'start', behavior: 'smooth' }); }
      else $('#pmap-card').scrollIntoView({ block: 'start', behavior: 'smooth' });
    }));
    $('#zoom-reset').addEventListener('click', () => { PV.zoom = null; drawProtMap(ctx); });
    $$('[data-lab]').forEach((b) => b.addEventListener('click', () => { PV.labels = b.dataset.lab; $$('[data-lab]').forEach((x) => x.setAttribute('aria-pressed', String(x === b))); drawProtMap(ctx); }));
    $$('[data-kind]').forEach((cb) => cb.addEventListener('change', () => { const s = mapKinds(); if (cb.checked) s.add(cb.dataset.kind); else s.delete(cb.dataset.kind); store.set('map-kinds', JSON.stringify([...s])); drawProtMap(ctx); }));
    $$('[data-row]').forEach((cb) => cb.addEventListener('change', () => { const s = mapRows(); if (cb.checked) s.add(cb.dataset.row); else s.delete(cb.dataset.row); store.set('map-rows', JSON.stringify([...s])); drawProtMap(ctx); }));
    const fig = (kind) => C.exportFigure(kind, protMapSvg(ctx, { exp: true }), geneOf(p) + '_idr_map' + (PV.zoom ? '_' + PV.zoom[0] + '-' + PV.zoom[1] : ''));
    $('#fig-svg').addEventListener('click', () => fig('svg'));
    $('#fig-png').addEventListener('click', () => fig('png'));
    $('#copy-fasta').addEventListener('click', (e) => { try { navigator.clipboard.writeText('>' + acc + ' ' + geneOf(p) + '\n' + seq.replace(/(.{60})/g, '$1\n') + '\n'); e.currentTarget.textContent = 'Copied'; } catch (err) { /* no clipboard */ } });
    $$('details.fold[data-fold^="p-"]').forEach((d) => d.addEventListener('toggle', () => UI.pOpen.open.set(d.dataset.fold, d.open)));
    // live: UniProt's function and regions, AlphaFold's confidence; the map redraws as each arrives
    C.loadUniProt(acc).then((u) => {
      if (token !== routeToken || !$('#pfunc')) return;
      ctx.uni = u;
      $('#pfunc').innerHTML = u && u.func ? '<span class="pf-text">' + esc(u.func.replace(/\s*\(PubMed:[^)]*\)/g, '')) + '</span> <span class="muted">Function, from <a href="https://www.uniprot.org/uniprotkb/' + acc + '/entry#function" target="_blank" rel="noopener">UniProt</a></span>' : '<span class="muted">' + (u ? 'No function text in UniProt.' : 'UniProt could not be reached.') + '</span>';
      drawProtMap(ctx);
    });
    let afMode = 'charge';
    const afFold = $('details[data-fold="p-3d"]');
    C.loadAF(acc).then((af) => { if (token !== routeToken) return; ctx.af = af; ctx.afTried = true; drawProtMap(ctx); if (afFold && afFold.open) structureView(ctx, afMode); });
    afFold.addEventListener('toggle', () => { if (afFold.open) structureView(ctx, afMode); });
    $$('[data-afc]').forEach((b) => b.addEventListener('click', () => { afMode = b.dataset.afc; $$('[data-afc]').forEach((x) => x.setAttribute('aria-pressed', String(x === b))); structureView(ctx, afMode); }));
    // the folds' tables: a row shows its place on the map
    const goSeg = (a, b, what, idr) => { PV.zoom = n > 300 ? [Math.max(1, a - 15), Math.min(n, b + 15)] : null; PV.hl = [a, b]; PV.pick = { seg: [a, b], what, idr }; drawProtMap(ctx); pickPanel(ctx); $('#pmap-card').scrollIntoView({ block: 'start', behavior: 'smooth' }); };
    const fl = Rp.fl || [], stF = { key: 'fs', dir: 1, n: 200 }, stS = { key: 'pos', dir: 1, n: 100 };
    const drawFl = () => tbl($('#pfl'), { id: 'p-fl', caption: 'Domain flanks of ' + geneOf(p), rows: fl.map((_, k) => k), st: stF, onState: drawFl, note: 'Every domain edge’s flank; click one to see it on the map.',
      onOpen: (k) => goSeg(fl[k][1], fl[k][2], (fl[k][0] ? 'C' : 'N') + '-side flank'),
      cols: [{ k: 'fs', label: 'Flank', tdCls: 'nw', get: (k) => fl[k][1], td: (k) => fl[k][1] + '–' + fl[k][2] + ' <span class="muted">' + (fl[k][0] ? 'after' : 'before') + '</span>' },
        { k: 'dom', label: 'Domain', get: (k) => fl[k][3], td: (k) => fl[k][3] + '–' + fl[k][4] },
        { k: 'za', label: 'Acidic score', num: true, get: (k) => (fl[k][9] == null ? null : fl[k][9] / 10), td: (k) => pctTop(fl[k][9]) + TOPB((fl[k][10] || 0) & 1) },
        { k: 'zb', label: 'Basic score', num: true, get: (k) => (fl[k][8] == null ? null : fl[k][8] / 10), td: (k) => pctTop(fl[k][8]) + TOPB((fl[k][10] || 0) & 2) },
        { k: 'kap', label: 'Charge segregation (κ)', num: true, get: (k) => (fl[k][11] == null ? null : fl[k][11] / 1000), td: (k) => f3(fl[k][11]) },
        { k: 'chem', label: 'Chemistry', get: (k) => CHEM_SHORT[fl[k][7]] || null }] });
    const all = [];
    for (const [k, l] of MAP_KINDS) for (const v of S[k] || []) all.push({ k, l, v, pos: v[0] });
    const ck2s = new Set((S.ck2 || []).map((v) => v[0]));
    const segAt = (pos) => idrs.find((i) => ID.s[i] <= pos && ID.e[i] >= pos);
    const drawS = () => tbl($('#psites'), { id: 'p-sites', caption: 'Sites of ' + geneOf(p), rows: all.map((_, j) => j), st: stS, onState: drawS, note: 'Click a site to read it on the map.',
      onOpen: (j) => { PV.pick = { pos: all[j].pos }; PV.zoom = n > 300 ? [Math.max(1, all[j].pos - 60), Math.min(n, all[j].pos + 60)] : PV.zoom; drawProtMap(ctx); pickPanel(ctx); $('#pmap-card').scrollIntoView({ block: 'start', behavior: 'smooth' }); },
      cols: [{ k: 'pos', label: 'Site', num: true, get: (j) => all[j].pos, td: (j) => '<b>' + residueName(seq, all[j].pos) + '</b>' },
        { k: 'kind', label: 'Modification', get: (j) => all[j].l, td: (j) => '<i class="kdot" style="background:var(--' + KCOL[all[j].k === 'ph' && ck2s.has(all[j].pos) ? 'ck2' : all[j].k] + ')"></i>' + esc(all[j].l) + (all[j].k === 'ph' && ck2s.has(all[j].pos) ? ' <span class="tag">CK2</span>' : '') },
        { k: 'idr', label: 'Its IDR', get: (j) => { const i = segAt(all[j].pos); return i == null ? null : ID.ct[i]; }, td: (j) => { const i = segAt(all[j].pos); return i == null ? '<span class="muted">ordered</span>' : ctChip(ID.ct[i]) + ' <span class="muted">' + ID.s[i] + '–' + ID.e[i] + '</span>'; } },
        { k: 'dq', label: 'Charge change', num: true, get: (j) => { const q = modQ(all[j].k, seq[all[j].pos - 1]); return q == null ? null : q - (RES_Q[seq[all[j].pos - 1]] || 0); }, td: (j) => { const q = modQ(all[j].k, seq[all[j].pos - 1]); return q == null ? '<span class="muted">—</span>' : fsign(q - (RES_Q[seq[all[j].pos - 1]] || 0), 0); }, tip: 'The charge the modification adds to the residue (idrkit.ptm)' },
        { k: 'src', label: 'Reported by', get: (j) => srcText(all[j].v), td: (j) => srcLinks(all[j].v[1], all[j].v[2], acc) },
        { k: 'note', label: 'Detail', get: (j) => { const v = all[j].v; return all[j].k === 'ph' && v.length > 4 ? 'CK2 percentile ' + v[3] + (v[7] != null ? '; occupancy ' + (v[7] / 10).toFixed(1) + '%' : '') : all[j].k === 'kac' && typeof v[3] === 'number' ? 'stoichiometry ' + (100 * v[3]).toFixed(3) + '%' : typeof v[3] === 'string' ? v[3] : null; } }] });
    drawFl(); drawS();
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
    if (e.key === 'Enter') {
      await loadProteins();
      // "FUS S42" or "P35637 R216": the protein, at that site
      const site = input.value.trim().match(/^(\S+)\s+[A-Za-z]?(\d+)$/);
      const m = matches(site ? site[1] : input.value);
      if (m.length) location.hash = PR.acc[m[0]] + (site ? '?site=' + site[2] : '');
      sug.hidden = true;
    }
    if (e.key === 'ArrowDown' && !sug.hidden) { e.preventDefault(); const a = sug.querySelector('a'); if (a) a.focus(); }
    if (e.key === 'Escape') sug.hidden = true;
  });
  sug.addEventListener('click', () => { sug.hidden = true; input.value = ''; });
  document.addEventListener('click', (e) => { if (!e.target.closest('.search')) sug.hidden = true; });
  // the explorer's colour tokens (--t-*, its TOKENS), which its stylesheet reads: set as it sets them, with the theme
  C.applyTokens();
  // figures drawn with the theme's colours (the protein map, the charts) are drawn again when it changes
  try { matchMedia('(prefers-color-scheme: dark)').addEventListener('change', () => { C.applyTokens(); route(); }); } catch (e) { /* old browsers */ }
  $('#theme-btn').addEventListener('click', () => {
    const cur = document.documentElement.getAttribute('data-theme') || (matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light');
    const next = cur === 'dark' ? 'light' : 'dark';
    document.documentElement.setAttribute('data-theme', next);
    try { localStorage.setItem('idr-theme', next); } catch (e) { /* private mode */ }
    C.applyTokens(next);
    route();
  });
  loadBase().then(() => { $('#built').textContent = 'Data built ' + B.generated + '.'; }).catch(() => {});
  route();
})();
