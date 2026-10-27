/* IDR Explorer. Data: docs/data/, written by Analysis/scripts/build_idr_explorer.py. The per-protein files use the
   Methylarginine Site Explorer's 64 shards and its hash (PTM_Pipeline docs/assets/explorer.js shardOf and
   scripts/36c_export_interpro_domains.shard), and its domain rows; its charge track definition (EMBOSS, 5-residue
   mean) is the one drawn here. */
(() => {
  'use strict';

  const $ = (s, el = document) => el.querySelector(s);
  const esc = (s) => String(s ?? '').replace(/[&<>"']/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
  const num = (x) => (x == null ? '' : Number(x).toLocaleString('en-US'));
  const SHARDS = 64;
  const shardOf = (acc) => { let h = 0; for (const c of acc) h = (h * 31 + c.charCodeAt(0)) % 1000003; return h % SHARDS; };
  const PTM_SITE = 'https://shechterlab.github.io/PTM_Pipeline/';
  const GH = {
    ptm: 'https://github.com/Shechterlab/PTM_Pipeline/blob/',
    site: 'https://github.com/Shechterlab/idr-explorer/blob/main/',
  };
  const PAGE = 100;

  // ---------------------------------------------------------------- data
  const cache = new Map();
  const getJSON = (path) => {
    if (!cache.has(path)) {
      cache.set(path, fetch(path).then((r) => { if (!r.ok) throw new Error(path + ': HTTP ' + r.status); return r.json(); }));
    }
    return cache.get(path);
  };
  let B = null, PR = null, FL = null, ID = null;
  const byAcc = new Map(), byGene = new Map();
  const loadBuild = async () => (B = B || await getJSON('data/build.json'));
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

  const z = (v) => (v == null ? null : v / 100);
  const fz = (v) => (v == null ? '' : (v >= 0 ? '+' : '−') + Math.abs(v / 100).toFixed(2));
  const zChip = (v, pole) => {
    if (v == null) return '<span class="muted">–</span>';
    const a = Math.max(0, Math.min(1, (v / 100) / 3));
    return '<span class="z" style="background:color-mix(in srgb, var(--' + pole + ') ' + Math.round(a * 55) + '%, transparent)">' + fz(v) + '</span>';
  };
  const CHEM_SHORT = ['Unmarked', 'R-type', 'K-type', 'Hybrid'];
  const FLAGS = [[1, 'TF'], [2, 'RBD'], [8, 'chaperone'], [16, 'HMG'], [4, 'histone-like']];
  const flagText = (f) => FLAGS.filter(([b]) => f & b).map(([, t]) => t).join(', ');
  const geneOf = (p) => PR.gene[p] || PR.acc[p];
  const protLink = (p, extra) => '<a href="#' + esc(PR.acc[p]) + (extra ? '?' + extra : '') + '">' + esc(geneOf(p)) + '</a>';

  // ---------------------------------------------------------------- routing
  const parseHash = () => {
    let h = location.hash.slice(1);
    try { h = decodeURIComponent(h); } catch (e) { /* keep as given */ }
    const [path, qs] = h.split('?');
    return { path: path || 'ideas', params: new URLSearchParams(qs || '') };
  };
  const setParams = (path, params) => {
    const s = params.toString();
    history.replaceState(null, '', '#' + path + (s ? '?' + s : ''));
  };
  const view = () => $('#view');
  const setTab = (name) => document.querySelectorAll('.tabs a').forEach((a) => {
    if (a.dataset.view === name) a.setAttribute('aria-current', 'page'); else a.removeAttribute('aria-current');
  });
  let routeToken = 0;
  async function route() {
    const token = ++routeToken;
    const { path, params } = parseHash();
    const v = view();
    const views = { ideas: ideasView, flanks: flanksView, idrs: idrsView, about: aboutView };
    setTab(views[path] ? path : '');
    if (views[path]) document.title = 'IDR Explorer · Shechter Lab';
    v.innerHTML = '<p class="loading">Loading…</p>';
    try {
      await loadBuild();
      if (token !== routeToken) return;
      if (views[path]) await views[path](params, token);
      else await proteinView(path, params, token);
    } catch (e) {
      if (token === routeToken) v.innerHTML = '<div class="panel"><b>Could not load this view.</b> <span class="muted">' + esc(e.message) + '</span></div>';
      console.error(e);
    }
  }
  window.addEventListener('hashchange', route);

  // ---------------------------------------------------------------- shared table
  // A sortable, paginated table over an array of row indices; columns: {key, label, num, get(i), html(i), sort(i)}.
  function table(host, opts) {
    const { caption, columns, rows, sortKey, desc, rowHeader = 0, csvName } = opts;
    let { page } = opts;
    // without onState the table keeps its own sort and page
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
      return '<th scope="col"' + (sorted ? ' aria-sort="' + sorted + '"' : '') + '><button type="button" data-sort="' + c.key + '">' + esc(c.label) + '</button></th>';
    }).join('');
    const body = slice.map((i) => '<tr>' + columns.map((c, k) => {
      const cell = c.html ? c.html(i) : esc(c.get(i) ?? '');
      return k === rowHeader ? '<th scope="row" data-label="' + esc(c.label) + '">' + cell + '</th>'
        : '<td' + (c.num ? ' class="num"' : '') + ' data-label="' + esc(c.label) + '">' + cell + '</td>';
    }).join('') + '</tr>').join('');
    host.innerHTML =
      '<div class="table-wrapper" tabindex="0" role="region" aria-labelledby="' + id + '"><table class="data"><caption id="' + id + '">' + esc(caption) +
      '</caption><thead><tr>' + head + '</tr></thead><tbody>' + body + '</tbody></table></div>' +
      '<div class="pager"><button class="btn" type="button" data-page="' + (page - 1) + '"' + (page <= 1 ? ' disabled' : '') + '>Previous</button>' +
      '<span>Page ' + page + ' of ' + num(pages) + ' · ' + num(order.length) + ' rows</span>' +
      '<button class="btn" type="button" data-page="' + (page + 1) + '"' + (page >= pages ? ' disabled' : '') + '>Next</button>' +
      (csvName ? '<button class="btn" type="button" data-csv="1">Download these ' + num(order.length) + ' rows (TSV)</button>' : '') + '</div>';
    host.querySelectorAll('[data-sort]').forEach((b) => b.addEventListener('click', () => {
      const k = b.dataset.sort;
      onState({ sort: k, desc: k === col.key ? !desc : true, page: 1 });
    }));
    host.querySelectorAll('[data-page]').forEach((b) => b.addEventListener('click', () => onState({ page: +b.dataset.page })));
    const csv = host.querySelector('[data-csv]');
    if (csv) {
      csv.addEventListener('click', () => {
        const lines = [columns.map((c) => c.label).join('\t')].concat(order.map((i) => columns.map((c) => String(c.get(i) ?? '')).join('\t')));
        const url = URL.createObjectURL(new Blob([lines.join('\n') + '\n'], { type: 'text/tab-separated-values' }));
        const a = document.createElement('a');
        a.href = url; a.download = csvName; document.body.append(a); a.click(); a.remove();
        setTimeout(() => URL.revokeObjectURL(url), 1000);
      });
    }
  }

  // A binned density of two columns (all rows in grey, the selection in the accent), as SVG.
  function density({ xs, ys, sel, xLabel, yLabel, xr, yr, bins = 56, marks = [] }) {
    const W = 420, H = 340, L = 46, R = 10, T = 10, Bm = 40;
    const pw = W - L - R, ph = H - T - Bm;
    const bx = (x) => Math.min(bins - 1, Math.max(0, Math.floor((x - xr[0]) / (xr[1] - xr[0]) * bins)));
    const by = (y) => Math.min(bins - 1, Math.max(0, Math.floor((y - yr[0]) / (yr[1] - yr[0]) * bins)));
    const grid = (idx) => {
      const g = new Uint32Array(bins * bins);
      for (const i of idx) { const x = xs[i], y = ys[i]; if (x == null || y == null) continue; g[by(y) * bins + bx(x)]++; }
      return g;
    };
    const all = grid(xs.keys());
    const rect = (g, color, scale) => {
      let max = 0; for (const v of g) if (v > max) max = v;
      let out = '';
      const cw = pw / bins, ch = ph / bins;
      for (let k = 0; k < g.length; k++) {
        if (!g[k]) continue;
        const a = 0.12 + 0.88 * Math.log1p(g[k]) / Math.log1p(max || 1);
        const gx = k % bins, gy = Math.floor(k / bins);
        out += '<rect x="' + (L + gx * cw).toFixed(1) + '" y="' + (T + ph - (gy + 1) * ch).toFixed(1) + '" width="' + (cw + 0.3).toFixed(1) + '" height="' + (ch + 0.3).toFixed(1) +
          '" style="fill:' + color + ';fill-opacity:' + (a * scale).toFixed(2) + '"><title>' + g[k] + '</title></rect>';
      }
      return out;
    };
    const sx = (x) => L + (x - xr[0]) / (xr[1] - xr[0]) * pw;
    const sy = (y) => T + ph - (y - yr[0]) / (yr[1] - yr[0]) * ph;
    const ticks = (r) => { const out = []; for (let t = Math.ceil(r[0]); t <= Math.floor(r[1]); t++) out.push(t); return out; };
    let svg = '<svg class="plane" viewBox="0 0 ' + W + ' ' + H + '" role="img" aria-label="' + esc(yLabel + ' against ' + xLabel) + '">';
    svg += '<rect x="' + L + '" y="' + T + '" width="' + pw + '" height="' + ph + '" style="fill:var(--surface-2)"/>';
    svg += rect(all, 'var(--muted)', sel ? 0.55 : 1);
    if (sel) svg += rect(grid(sel), 'var(--accent)', 1);
    for (const t of ticks(xr)) svg += '<line x1="' + sx(t) + '" x2="' + sx(t) + '" y1="' + (T + ph) + '" y2="' + (T + ph + 4) + '" style="stroke:var(--muted)"/><text x="' + sx(t) + '" y="' + (T + ph + 16) + '" text-anchor="middle" font-size="12">' + t + '</text>';
    for (const t of ticks(yr)) svg += '<line x1="' + (L - 4) + '" x2="' + L + '" y1="' + sy(t) + '" y2="' + sy(t) + '" style="stroke:var(--muted)"/><text x="' + (L - 7) + '" y="' + (sy(t) + 4) + '" text-anchor="end" font-size="12">' + t + '</text>';
    if (xr[0] < 0 && xr[1] > 0) svg += '<line x1="' + sx(0) + '" x2="' + sx(0) + '" y1="' + T + '" y2="' + (T + ph) + '" style="stroke:var(--line-strong);stroke-dasharray:3 3"/>';
    if (yr[0] < 0 && yr[1] > 0) svg += '<line x1="' + L + '" x2="' + (L + pw) + '" y1="' + sy(0) + '" y2="' + sy(0) + '" style="stroke:var(--line-strong);stroke-dasharray:3 3"/>';
    for (const m of marks) {
      svg += '<circle cx="' + sx(m.x).toFixed(1) + '" cy="' + sy(m.y).toFixed(1) + '" r="3.5" style="fill:none;stroke:var(--ink);stroke-width:1.4"/>' +
        '<text x="' + (sx(m.x) + 5).toFixed(1) + '" y="' + (sy(m.y) - 4).toFixed(1) + '" font-size="11" style="fill:var(--ink)">' + esc(m.label) + '</text>';
    }
    svg += '<text x="' + (L + pw / 2) + '" y="' + (H - 6) + '" text-anchor="middle" font-size="12">' + esc(xLabel) + '</text>';
    svg += '<text transform="translate(14 ' + (T + ph / 2) + ') rotate(-90)" text-anchor="middle" font-size="12">' + esc(yLabel) + '</text>';
    return svg + '</svg>';
  }

  // ---------------------------------------------------------------- ideas
  const src = (repo, path, label) => '<a href="' + GH[repo] + path + '">' + esc(label || path.split('/').pop()) + '</a>';
  // MicroDomains and glu_proteomics are private until publication: their files are named, not linked
  const mdsrc = (path) => '<span class="src" title="Shechterlab/MicroDomains (private until publication)">MicroDomains: ' + esc(path) + '</span>';
  const glusrc = (path) => '<span class="src" title="Shechterlab/glu_proteomics (private until publication)">glu_proteomics: ' + esc(path) + '</span>';

  function anatomy() {
    // a folded domain, its two edges, and a dtIDR of each pole
    return '<svg class="anatomy" viewBox="0 0 760 150" role="img" aria-label="A folded domain with a basic dtIDR on one edge and an acidic dtIDR on the other">' +
      '<line x1="10" y1="75" x2="750" y2="75" style="stroke:var(--line-strong);stroke-width:3"/>' +
      '<rect x="300" y="45" width="160" height="60" rx="12" style="fill:var(--domain);stroke:var(--domain-ink)"/>' +
      '<text x="380" y="80" text-anchor="middle" font-size="14" style="fill:var(--domain-ink)">folded domain</text>' +
      '<rect x="150" y="66" width="150" height="18" rx="4" style="fill:var(--basic-wash);stroke:var(--basic)"/>' +
      '<text x="225" y="79" text-anchor="middle" font-size="12" class="mono" style="fill:var(--basic)">…RGRGRK…</text>' +
      '<rect x="460" y="66" width="150" height="18" rx="4" style="fill:var(--acidic-wash);stroke:var(--acidic)"/>' +
      '<text x="535" y="79" text-anchor="middle" font-size="12" class="mono" style="fill:var(--acidic)">…EEDEEDEE…</text>' +
      '<text x="225" y="52" text-anchor="middle" font-size="12">basic dtIDR (≤50 aa)</text>' +
      '<text x="225" y="104" text-anchor="middle" font-size="11" class="lab" style="fill:var(--muted)">PRMT / KMT / KAT writers</text>' +
      '<text x="535" y="52" text-anchor="middle" font-size="12">acidic dtIDR (≤50 aa)</text>' +
      '<text x="535" y="104" text-anchor="middle" font-size="11" style="fill:var(--muted)">CK2, TTLL writers</text>' +
      '<path d="M260 118 h40" style="stroke:var(--muted)"/><text x="280" y="134" text-anchor="middle" font-size="10" style="fill:var(--muted)">~25 aa: few marks</text>' +
      '</svg>';
  }

  async function ideasView() {
    const c = B.counts, h = B.highlights, X = B.ideas;
    const pct = (x) => (100 * x).toFixed(x > 0.999 ? 2 : 1);
    const p = (v) => (v < 0.001 ? 'p < 0.001' : 'p = ' + (v < 0.01 ? v.toFixed(3) : v.toFixed(2)));
    const or = (v) => v.toFixed(2);
    const fig = (name, alt, cap) => '<figure class="fig"><a href="figures/' + name + '.png"><img src="figures/' + name + '.png" alt="' + esc(alt) +
      '" loading="lazy"></a>' + (cap ? '<figcaption>' + cap + '</figcaption>' : '') + '</figure>';
    const pairRow = (r, a, b) => '<tr><td><a href="#' + r[1] + '?f=' + r[3] + '-' + r[4] + '">' + esc(r[0]) + '</a></td><td>' + r[2] + '-side ' + r[3] + '–' + r[4] +
      '</td>' + (a ? '<td>' + esc(r[5]) + '</td>' : '') + '<td class="num">' + r[6].toFixed(2) + '</td>' + (b ? '<td class="num">' + pct(r[7]) + '</td>' : '') + '</tr>';
    const A = X.acidic, F = X.families, T = X.trans, PT = X.patterning.tests, PG = X.patterning.groups, CI = X.cis;
    const net = T.networks.string700;
    const chem = h.chemistry_counts;
    const ideas = [
      {
        cls: 'lead', tag: ['data', 'The unit'], title: 'A dtIDR is an address, not a property',
        body: '<p>“IDRs carry modifications” is a statement about a category of sequence. A domain-tethered IDR is a bounded object: one folded domain, one side, and the residues beside it, at coordinates someone can clone. Every domain edge of the reviewed human proteome is scored from sequence alone — disorder, net charge, local R/K clustering and RG/GR density — and no modification enters the score. Read the same windows with the charge inverted and you get the second pole.</p>' + anatomy() +
          '<ul><li><b>' + num(c.flank_windows) + '</b> distinct windows on <b>' + num(c.proteins_with_flanks) + '</b> proteins, and <b>' + num(c.idrs) + '</b> IDRs of 20 or more residues on ' + num(c.proteins_with_idr) + ' proteins.</li>' +
          '<li>The basic score finds the windows that carry an annotated methylarginine: AUC 0.701, odds ratio 2.13 per SD (2.02–2.23) over the full proteome — the same as in the earlier PTM-selected universe.</li>' +
          '<li>Disorder alone does worst (AUC 0.601); net charge 0.649, RG/GR density 0.670.</li></ul>',
        foot: [mdsrc('Analysis/microdomains/full_proteome_migration.md'), '<a href="#flanks?sort=zb">Rank the flanks →</a>', '<a href="#idrs">Every IDR →</a>'],
      },
      {
        tag: ['data', 'Position'], title: 'Marks avoid the first ~25 residues beside the fold',
        body: '<p>Methyl-R, methyl-K and acetyl-K sites outside domains are depleted within 10 residues of a domain edge (0.80× of a within-protein null, empirical p = 0.004) and reach baseline about 30 residues out. The functional E2F1 and BRD4 sites sit inside that zone, so the claim is about consequence, not density: a mark there can reach the neighbouring fold.</p>',
        foot: [mdsrc('Analysis/microdomains/why_dtidrs_are_not_just_idrs.md'), '<a href="#E2F1">E2F1 →</a>', '<a href="#BRD4">BRD4 →</a>'],
      },
      {
        tag: ['data', 'Classes'], title: 'The classes that need a tuner have one',
        body: '<table class="mini"><tr><th>Domain-edge class</th><th class="num">mean z</th></tr>' +
          [['Transcription factor / DNA binding', '+0.444'], ['RNA binding', '+0.365'], ['Chromatin reader', '+0.287'], ['Signalling & catalytic', '+0.056'], ['Protein-interaction scaffold', '−0.006'], ['Membrane & extracellular', '−0.292']]
            .map(([x, y]) => '<tr><td>' + x + '</td><td class="num">' + y + '</td></tr>').join('') + '</table>' +
          '<p class="caveat">Disorder alone separates these classes about as well; what the composite adds is where the methylation is.</p>',
        foot: [mdsrc('Analysis/microdomains/full_proteome_migration.md'), '<a href="#flanks?sort=zb&tf=1">TF flanks →</a>'],
      },
      {
        tag: ['data', 'Chemistry'], title: 'Two basic chemistries, rarely on one flank',
        body: '<p>R-type flanks carry a GAR/RGG array (PRMT grammar); K-type flanks a lysine mark (SETD6, p300). With annotation depth and the score in the model, a GAR array predicts <i>against</i> a lysine mark: odds ratio 0.26 (0.16–0.44) over 120,870 flanks, 0.43 (0.30–0.62) among studied proteins. Hybrids are rare — ' + num(chem['Hybrid (both chemistries)']) +
          ' flanks against ' + num(chem['R-type (GAR array)']) + ' R-type and ' + num(chem['K-type (lysine-marked)']) + ' K-type — and E2F1 and BRD4 are both hybrids.</p>' +
          '<table class="mini"><tr><th>Protein</th><th>Best basic flank</th><th class="num">z</th><th class="num">percentile</th></tr>' + h.focal_basic.map((r) => pairRow(r, false, true)).join('') + '</table>',
        foot: [mdsrc('Analysis/microdomains_review_bundle/tables/microdomain_chemistry_exclusivity.tsv'), '<a href="#flanks?chem=3&sort=zb">The hybrids →</a>'],
      },
      {
        cls: 'wide', tag: ['new', 'Acidic pole, full proteome'], title: 'The charge mirror: histone chaperones and CK2',
        body: '<p>Invert the charge term and replace RG/GR with glutamate-in-an-acidic-run, on the same ' + num(c.flank_windows) + ' windows. Top decile against the rest: CK2-attributed phosphorylation OR <b>' + or(A.ck2[0]) + '</b> (' + A.ck2[2] + ' flanks), and ' + or(A.ck2_share[0]) +
          ' among flanks that are phosphorylated at all, where the basic pole gives ' + or(A.basic_ck2_share[0]) + '; histone chaperones ' + or(A.chaperone[0]) + '; arginine methylation ' + or(A.methylR[0]) + ' (the basic pole: ' + or(A.basic_methylR[0]) +
          '). The NAP-like chaperone fold is the most acidic family, ' + or(F.nap.acidic_or_per_sd) + ' per SD (' + or(F.nap.acidic_ci_low) + '–' + or(F.nap.acidic_ci_high) + '), flat at the basic pole (' + or(F.nap.basic_or_per_sd) +
          '). Chaperones that the earlier, methylated-protein universe left out now rank near the top: ' + h.chaperones_full_proteome.map((r) => '<a href="#' + r[1] + '?f=' + r[3] + '-' + r[4] + '">' + esc(r[0]) + '</a> ' + pct(r[7])).join(', ') + ' (percentile).</p>' +
          fig('figH', 'The two poles: density, writers, families and chaperones', 'Two axes, not one; each pole has its own writer; the chaperone fold is acidic, the histone fold basic.') +
          '<p class="caveat">Did not survive the move to the full proteome: the histone fold is no longer depleted at the acidic pole (' + or(F.histone_fold.acidic_or_per_sd) + ' per SD, q = ' + F.histone_fold.acidic_q_value.toFixed(2) + '; it was 0.60) and histone-like proteins are flat (' + or(A.histone_like[0]) + '). The substrate is still basic: the histone fold is ' + or(F.histone_fold.basic_or_per_sd) + ' per SD at the basic pole.</p>' +
          '<table class="mini"><tr><th>Gene</th><th>Flank</th><th>Adjacent domain</th><th class="num">z</th></tr>' + h.top_acidic.slice(0, 8).map((r) => pairRow(r, true, false)).join('') + '</table>',
        foot: [mdsrc('Analysis/microdomains/acidic_dtidr_notes.md, section 0'), '<a href="#flanks?sort=za">Rank the acidic pole →</a>', '<a href="#flanks?sort=za&chap=1">Chaperones and HMG →</a>'],
      },
      {
        tag: ['new', 'Patterning'], title: 'Blocky at the acidic pole, mixed at the basic',
        body: '<p>Every IDR scored for local charge blocks (Das & Pappu’s δ, as a z-score against shuffles of its own residues). Against the other IDRs of the same protein, with length and charge fraction held: IDRs with a CK2 site are <b>' + fz(PT.ck2[0] * 100) + '</b> SD blockier (' + p(PT.ck2[1]) + ', ' + PT.ck2[2] + ' proteins), with any phosphosite ' + fz(PT.ph[0] * 100) +
          ', while methylated IDRs are better mixed (' + fz(PT.mR[0] * 100) + ', ' + p(PT.mR[1]) + ') and so are acyl-lysine ones (' + fz(PT.kac[0] * 100) + '). Glutamylated IDRs point the blocky way on ' + PT.glu[2] + ' proteins (' + fz(PT.glu[0] * 100) + ', ' + p(PT.glu[1]) +
          '). The positive control, DDX4 1–236, whose alternating charge blocks drive its phase separation, scores ' + PG['positive control: DDX4 1-236'][0].toFixed(1) + '.</p>' +
          fig('figJ', 'Methylarginines and acidic tracts in cis; charge patterning of every IDR; modified IDRs within protein', 'Right: change in δ z of IDRs carrying each mark against the other IDRs of the same protein.'),
        foot: [mdsrc('Analysis/microdomains/idr_charge_patterning_notes.md'), '<a href="#idrs?sort=dz">The blockiest IDRs →</a>'],
      },
      {
        tag: ['new', 'In cis'], title: 'The two charges share proteins but not segments',
        body: '<p>Methylated proteins carry a disordered acidic tract (a run of ≥6 D/E) more often than other reviewed proteins: ' + CI.protein6[0] + ' of ' + num(CI.protein6[1]) + ' against ' + num(CI.protein6[2]) + ' of ' + num(CI.protein6[3]) +
          ', Mantel–Haenszel OR ' + or(CI.protein6[4]) + ' (' + or(CI.protein6[5]) + '–' + or(CI.protein6[6]) + ') within abundance deciles. Yet within those proteins the methylarginines keep away from the tracts: OR ' + or(CI.cis6_30[0]) + ' (' + or(CI.cis6_30[1]) + '–' + or(CI.cis6_30[2]) +
          ') for lying within 30 residues of one, and ' + or(CI.cis_range[0]) + '–' + or(CI.cis_range[1]) + ' over tract lengths 5–8 and distances 10–100 (rotation p ≤ ' + CI.cis_range[2].toFixed(2) + '). Without RG/GR arginines: ' + or(CI.without_rg_range[0]) + '–' + or(CI.without_rg_range[1]) + '.</p>',
        foot: [src('ptm', 'ccr-b90e3503-xnxh85/scripts/52_preview_acidic_tracts.py', 'PTM_Pipeline scripts/52'), '<a href="#idrs?mixed=1&sort=dz">IDRs with both charges →</a>'],
      },
      {
        tag: ['neg', 'Tested: not a general rule'], title: 'The poles do not pair across binding partners',
        body: '<p>If acidic dtIDRs read basic microdomains in trans, acidic–basic protein pairs should be over-represented among binding partners. Against classes shuffled within degree and nuclear status (STRING physical ≥ 700, histones left out), they are not: acidic–basic O/E ' + or(net['acidic-basic'][0]) + ', while basic–basic is ' + or(net['basic-basic'][0]) +
          ' — like binds like (complementarity ' + or(net['contrast: acidic-basic vs same-charge'][0]) + ', ' + p(net['contrast: acidic-basic vs same-charge'][1]) + '; the same in Complex Portal). Within proteins, IntAct regions on an acidic window bind partners only ' + T.regions.on_acidic_window.delta_partner_basic.toFixed(2) + ' SD more basic (p = ' + T.regions.on_acidic_window.p_partner_basic.toFixed(2) + ').</p>' +
          fig('figI', 'Pole pairing across STRING, Complex Portal and IntAct binding regions', 'Charge pairing is real in named systems such as chaperone and histone; it is not who-binds-whom across the proteome.'),
        foot: [mdsrc('Analysis/microdomains/acidic_dtidr_notes.md, section 10')],
      },
      {
        tag: ['data', 'Glutamylation'], title: 'Glutamylation follows acidity, not RNA binding',
        body: '<p>Re-analysed TTLL4 SILAC pulldowns (PRIDE PXD065251): E and D are the only enriched residues, the longest E run predicts enrichment (OR 1.87 per SD) and RNA binding adds nothing once acidity and abundance are in the model. The acidic dtIDR score predicts enrichment in every TTLL pulldown (TTLL4 OR 1.93 per SD; the TTLL8 glycylase 2.84), the basic score in none.</p>',
        foot: [glusrc('docs/ANALYSIS_NOTES.md, sections 0.3 and 0.10'), '<a href="#idrs?glu=1">IDRs with Glu sites →</a>'],
      },
      {
        tag: ['data', 'Where the poles meet'], title: 'Glutamate beside a basic patch',
        body: '<p>The TTLL4-dependent sites seen directly in the spectra sit beside basic patches: 14 of 15 have a patch (net ≥ +3 in 10 residues) within 30 residues, against 57% of covered glutamates. HNRNPU E262 follows <span class="mono">KKRGVKRPR</span> in a flank high on both poles; the YBX1 and YBX3 sites lie in alternating charge blocks. Bona fide sites: ' +
          ['HNRNPU', 'NCL', 'YBX3'].map((g) => '<a href="#' + g + '">' + g + '</a>').join(', ') + '; ' + num(c.glu_sites) + ' glutamylation sites shown in all, with NPM1 E126 and UniProt’s tubulin and CGAS annotations.</p>' +
          '<p class="caveat">A glutamylated peptide ionizes better beside K/R, so detection favours these sites; long acidic runs are nearly invisible to tryptic MS.</p>',
        foot: [glusrc('docs/figures/18_glutamylated__ttll4_site_list.tsv')],
      },
      {
        tag: ['hyp', 'Hypothesis'], title: 'A dtIDR is an avidity tuner',
        body: '<p>The folded domain chooses the site; the adjacent ~30 residues make a second, low-specificity contact with the same ligand surface and set how long the domain stays. Modification retunes that contact, so it should change residence time rather than site selection. The trans test sharpens it: a basic dtIDR’s partner is the polyanion it slides along (DNA, RNA, an acidic tail in cis), not an acidic dtIDR on another protein.</p>',
        foot: [mdsrc('Analysis/microdomains/hypothesis_and_next_experiments.md')],
      },
      {
        tag: ['new', 'Public data'], title: 'What the site integrates, and what is next',
        body: '<ul><li><b>UniProt 2026_03</b> modified residues (CC BY 4.0): ' + num(c.ck2_phosphosites) + ' CK2-attributed phosphosites with iPTMnet’s, ' + num(c.adp_ribosyl) + ' ADP-ribosylations, ' + num(c.citrulline) + ' citrullines, tubulin and CGAS polyglutamylation.</li>' +
          '<li><b>PTMeXchange</b> reprocessed sites (CC0), graded Gold/Silver/Bronze: acetyl, ubiquitin (' + num(c.ubiquitin_lysines) + ' lysines with iPTMnet’s) and SUMO (' + num(c.sumo_lysines) + '). Of its ' + num(h.ptmexchange_methyl_r.sites) + ' methylarginines, ' + num(h.ptmexchange_methyl_r.in_release) + ' are already in the methylarginine release.</li>' +
          '<li>Next, by value: kinase-motif scores for every S/T from the Ser/Thr kinome atlas (Johnson et al. 2023; CK2 among 303 kinases); a 2024 citrullinome (Rebak et al., 14,056 sites, 186 on methylated arginines); ADP-ribosylomes (Hendriks et al. 2021); FINCHES predictions of IDR–IDR attraction.</li></ul>',
        foot: ['<a href="#about">Sources and licences →</a>'],
      },
      {
        tag: ['open', 'Open'], title: 'What is not settled',
        body: '<ul><li>Why basic-pole proteins bind each other more than chance, after degree, nuclear status and annotated RNA binding are held.</li>' +
          '<li>Whether CK2 writes the blocky acidic IDRs it is enriched in, or only marks them: a motif-score test on every S/T is the next step.</li>' +
          '<li>Most TTLL4 sites are region-level; there is no proteome-wide glutamylation site table.</li>' +
          '<li>The independent interface test of the dtIDR score is not yet predictive (held-out AUC 0.536).</li></ul>',
        foot: [mdsrc('Analysis/microdomains/README.md')],
      },
    ];
    view().innerHTML =
      '<p class="lede">Disordered regions of the reviewed human proteome, read two ways: as <b>IDRs</b> (every MobiDB disordered region of 20 or more residues) and as <b>dtIDRs</b>, the segments tethered to the edge of a folded domain, scored for a basic pole that arginine and lysine writers act on and an acidic pole that CK2 and the glutamylases act on.</p>' +
      '<div class="stats">' + [
        [c.idrs, 'IDRs (≥20 aa)'], [c.flank_windows, 'dtIDR windows'], [c.acidic_stretches, 'acidic stretches'], [c.methylarginines, 'methylarginines'],
        [c.phosphosites, 'phosphosites (' + num(c.ck2_phosphosites) + ' CK2)'], [c.acyl_lysines + c.ubiquitin_lysines + c.sumo_lysines, 'acyl, ubiquitin and SUMO lysines'], [c.glu_sites, 'glutamylation sites'],
      ].map(([n, l]) => '<div class="stat"><b>' + num(n) + '</b><span>' + esc(l) + '</span></div>').join('') + '</div>' +
      '<div class="ideas">' + ideas.map((d) =>
        '<article class="idea ' + (d.cls || '') + '"><span class="tag ' + d.tag[0] + '">' + esc(d.tag[1]) + '</span><h3>' + esc(d.title) + '</h3>' + d.body +
        '<div class="foot">' + d.foot.join('') + '</div></article>').join('') + '</div>';
  }

  // ---------------------------------------------------------------- dtIDR flanks
  async function flanksView(params) {
    view().innerHTML = '<p class="loading">Loading 121k flanks…</p>';
    await Promise.all([loadProteins(), loadFlanks()]);
    const n = FL.p.length;
    const st = {
      q: params.get('q') || '', sort: params.get('sort') || 'za', desc: params.get('asc') !== '1', page: +(params.get('page') || 1),
      chem: params.get('chem') ?? '', side: params.get('side') || '', tf: params.get('tf') === '1', rbd: params.get('rbd') === '1',
      chap: params.get('chap') === '1', hist: params.get('hist') === '1', mr: params.get('mr') === '1', glu: params.get('glu') === '1',
      ck: params.get('ph') === '1', ck2: params.get('ck2') === '1',
    };
    view().innerHTML =
      '<h2>dtIDR flanks</h2><p class="lede">Every window beside a folded domain, with the basic score (<code>dtidr_score_z</code>) and the acidic score recomputed on the same windows. Modifications are counted inside the window and never enter a score.</p>' +
      '<div class="panel"><form class="filters" id="ff">' +
      '<label>Gene or accession<input type="text" name="q" value="' + esc(st.q) + '" placeholder="e.g. NCL"></label>' +
      '<label>Sort by<select name="sort">' + [['za', 'acidic score'], ['zb', 'basic score'], ['both', 'both poles (lower of the two)'], ['lar', 'longest D/E run'], ['gene', 'gene']]
        .map(([v, l]) => '<option value="' + v + '"' + (st.sort === v ? ' selected' : '') + '>' + l + '</option>').join('') + '</select></label>' +
      '<label>Chemistry<select name="chem"><option value="">any</option>' + CHEM_SHORT.map((l, k) => '<option value="' + k + '"' + (st.chem === String(k) ? ' selected' : '') + '>' + l + '</option>').join('') + '</select></label>' +
      '<label>Side of domain<select name="side"><option value="">either</option><option value="N"' + (st.side === 'N' ? ' selected' : '') + '>N-side</option><option value="C"' + (st.side === 'C' ? ' selected' : '') + '>C-side</option></select></label>' +
      [['tf', 'TF flank'], ['rbd', 'RNA-binding flank'], ['chap', 'histone chaperone or HMG'], ['mr', 'methyl-R inside'], ['glu', 'Glu site inside'], ['ph', 'phosphosite inside'], ['ck2', 'CK2 site inside'], ['hist', 'include histone-like']]
        .map(([k, l]) => '<label class="check"><input type="checkbox" name="' + k + '"' + ((k === 'ph' ? st.ck : st[k]) ? ' checked' : '') + '>' + l + '</label>').join('') +
      '</form></div><div class="plot-row"><div id="ft"></div><div class="panel" id="fplot"></div></div>';
    const form = $('#ff');
    const q = st.q.trim().toUpperCase();
    const qAccs = q ? new Set(PR.acc.map((a, i) => ((PR.gene[i] || '').toUpperCase().startsWith(q) || a.startsWith(q) ? i : -1)).filter((i) => i >= 0)) : null;
    const rows = [];
    for (let i = 0; i < n; i++) {
      const f = FL.fl[i];
      if (!st.hist && (f & 4)) continue;
      if (qAccs && !qAccs.has(FL.p[i])) continue;
      if (st.chem !== '' && FL.chem[i] !== +st.chem) continue;
      if (st.side && FL.side[i] !== (st.side === 'N' ? 0 : 1)) continue;
      if (st.tf && !(f & 1)) continue;
      if (st.rbd && !(f & 2)) continue;
      if (st.chap && !(f & 24)) continue;
      if (st.mr && !FL.mR[i]) continue;
      if (st.glu && !FL.glu[i]) continue;
      if (st.ck && !FL.ph[i]) continue;
      if (st.ck2 && !FL.ck2[i]) continue;
      rows.push(i);
    }
    const base = st.hist ? n : FL.fl.reduce((k, f) => k + !(f & 4), 0);
    const filtered = rows.length < base;
    const columns = [
      { key: 'gene', label: 'Gene', get: (i) => geneOf(FL.p[i]), html: (i) => protLink(FL.p[i], 'f=' + FL.fs[i] + '-' + FL.fe[i]) },
      { key: 'acc', label: 'Accession', get: (i) => PR.acc[FL.p[i]] },
      { key: 'side', label: 'Side', get: (i) => (FL.side[i] ? 'C' : 'N') },
      { key: 'fs', label: 'Flank', get: (i) => FL.fs[i] + '-' + FL.fe[i], sort: (i) => FL.fs[i] },
      { key: 'dn', label: 'Adjacent domain', get: (i) => FL.domain_names[FL.dn[i]], html: (i) => esc((FL.domain_names[FL.dn[i]] || '').split(';')[0].slice(0, 48)) },
      { key: 'zb', label: 'Basic z', num: true, get: (i) => z(FL.zb[i]), html: (i) => zChip(FL.zb[i], 'basic') },
      { key: 'za', label: 'Acidic z', num: true, get: (i) => z(FL.za[i]), html: (i) => zChip(FL.za[i], 'acidic') },
      { key: 'both', label: 'Lower of the two', num: true, get: (i) => (FL.zb[i] == null || FL.za[i] == null ? null : Math.min(FL.zb[i], FL.za[i]) / 100), html: (i) => fz(Math.min(FL.zb[i], FL.za[i])) },
      { key: 'chem', label: 'Chemistry', get: (i) => CHEM_SHORT[FL.chem[i]] || '' },
      { key: 'lar', label: 'Longest D/E run', num: true, get: (i) => FL.lar[i] },
      { key: 'mR', label: 'methyl-R', num: true, get: (i) => FL.mR[i] },
      { key: 'kme', label: 'methyl-K', num: true, get: (i) => FL.kme[i] },
      { key: 'kac', label: 'acyl-K', num: true, get: (i) => FL.kac[i] },
      { key: 'ph', label: 'phospho', num: true, get: (i) => FL.ph[i] },
      { key: 'ck2', label: 'CK2', num: true, get: (i) => FL.ck2[i] },
      { key: 'glu', label: 'Glu', num: true, get: (i) => FL.glu[i] },
      { key: 'flags', label: 'Class', get: (i) => flagText(FL.fl[i]) },
    ];
    const save = (over) => {
      Object.assign(st, over);
      const p = new URLSearchParams();
      if (st.q) p.set('q', st.q);
      if (st.sort !== 'za') p.set('sort', st.sort);
      if (!st.desc) p.set('asc', '1');
      if (st.page > 1) p.set('page', st.page);
      if (st.chem !== '') p.set('chem', st.chem);
      if (st.side) p.set('side', st.side);
      for (const k of ['tf', 'rbd', 'chap', 'mr', 'glu', 'ck2', 'hist']) if (st[k]) p.set(k, '1');
      if (st.ck) p.set('ph', '1');
      setParams('flanks', p);
    };
    const draw = () => table($('#ft'), {
      caption: 'dtIDR flanks', columns, rows, sortKey: st.sort, desc: st.desc, page: st.page, csvName: 'dtidr_flanks.tsv',
      onState: (o) => { save(o); draw(); $('#ft').scrollIntoView({ block: 'nearest' }); },
    });
    draw();
    const xs = FL.zb.map(z), ys = FL.za.map(z);
    let xr = [Infinity, -Infinity], yr = [Infinity, -Infinity];
    for (let i = 0; i < n; i++) { if (xs[i] != null) { xr[0] = Math.min(xr[0], xs[i]); xr[1] = Math.max(xr[1], xs[i]); } if (ys[i] != null) { yr[0] = Math.min(yr[0], ys[i]); yr[1] = Math.max(yr[1], ys[i]); } }
    xr = [Math.floor(xr[0] * 2) / 2, Math.ceil(xr[1] * 2) / 2]; yr = [Math.floor(yr[0] * 2) / 2, Math.ceil(yr[1] * 2) / 2];
    const marks = (B.highlights.focal_basic.slice(0, 2).concat(newlyMarks())).map((r) => {
      const p = byAcc.get(r[1]);
      for (let i = 0; i < n; i++) if (FL.p[i] === p && FL.fs[i] === r[3] && FL.fe[i] === r[4]) return { x: xs[i], y: ys[i], label: r[0] };
      return null;
    }).filter(Boolean);
    $('#fplot').innerHTML = '<h3 style="margin:0 0 6px;font-size:1.05rem">The two poles</h3>' +
      density({ xs, ys, sel: filtered ? rows : null, xLabel: 'basic dtIDR score (z)', yLabel: 'acidic dtIDR score (z)', xr, yr, marks }) +
      '<p class="caveat">Binned density of all windows' + (filtered ? '; the current selection in orange' : '') + '. Scores live on a lattice (most windows are 50 aa), so a scatter would draw stripes.</p>';
    form.addEventListener('change', (e) => {
      const t = e.target;
      if (t.type === 'checkbox') save({ [t.name === 'ph' ? 'ck' : t.name]: t.checked, page: 1 });
      else if (t.name === 'sort') save({ sort: t.value, desc: t.value !== 'gene', page: 1 });
      else save({ [t.name]: t.value, page: 1 });
      route();
    });
    form.addEventListener('submit', (e) => e.preventDefault());
  }
  const newlyMarks = () => (B.highlights.chaperones_new_on_full_proteome || []).filter((r) => ['NPM2', 'DAXX', 'ASF1A'].includes(r[0]));

  // ---------------------------------------------------------------- IDRs
  async function idrsView(params) {
    view().innerHTML = '<p class="loading">Loading 33k IDRs…</p>';
    await Promise.all([loadProteins(), loadIdrs()]);
    const n = ID.p.length;
    const b = (k) => params.get(k) === '1';
    const st = {
      q: params.get('q') || '', sort: params.get('sort') || 'len', desc: params.get('asc') !== '1', page: +(params.get('page') || 1),
      minlen: +(params.get('minlen') || 0), term: params.get('term') || '', as: b('as'), rg: b('rg'), bp: b('bp'), mixed: b('mixed'),
      mr: b('mr'), glu: b('glu'), dt: b('dt'), ck2: b('ck2'), blocky: b('blocky'), ncpr: params.get('ncpr') || '',
    };
    view().innerHTML =
      '<h2>IDRs</h2><p class="lede">Every MobiDB disordered region of 20 or more residues, with its charge and charge patterning (δ z: local charge blocks against shuffles of its own residues; above 2 is blocky), its acidic stretches, RG/RGG clusters and basic patches (idrkit), its modification sites, and the best dtIDR window it overlaps.</p>' +
      '<div class="panel"><form class="filters" id="if">' +
      '<label>Gene or accession<input type="text" name="q" value="' + esc(st.q) + '" placeholder="e.g. NPM1"></label>' +
      '<label>Minimum length<input type="number" name="minlen" min="20" step="10" value="' + (st.minlen || '') + '" placeholder="20"></label>' +
      '<label>Net charge<select name="ncpr"><option value="">any</option><option value="neg"' + (st.ncpr === 'neg' ? ' selected' : '') + '>acidic (NCPR ≤ −0.1)</option><option value="pos"' + (st.ncpr === 'pos' ? ' selected' : '') + '>basic (NCPR ≥ +0.1)</option></select></label>' +
      '<label>Position<select name="term"><option value="">any</option>' + [['N', 'N-terminal'], ['C', 'C-terminal'], ['I', 'internal']].map(([v, l]) => '<option value="' + v + '"' + (st.term === v ? ' selected' : '') + '>' + l + '</option>').join('') + '</select></label>' +
      [['as', 'acidic stretch'], ['rg', 'RG/RGG cluster'], ['bp', 'basic patch'], ['mixed', 'acidic stretch and a basic patch or RG cluster'], ['blocky', 'blocky (δ z ≥ 2)'], ['mr', 'methyl-R'], ['ck2', 'CK2 site'], ['glu', 'Glu site'], ['dt', 'overlaps a dtIDR window']]
        .map(([k, l]) => '<label class="check"><input type="checkbox" name="' + k + '"' + (st[k] ? ' checked' : '') + '>' + l + '</label>').join('') +
      '</form></div><div class="plot-row"><div id="it"></div><div class="panel" id="iplot"></div></div>';
    const q = st.q.trim().toUpperCase();
    const qAccs = q ? new Set(PR.acc.map((a, i) => ((PR.gene[i] || '').toUpperCase().startsWith(q) || a.startsWith(q) ? i : -1)).filter((i) => i >= 0)) : null;
    const rows = [];
    for (let i = 0; i < n; i++) {
      const len = ID.e[i] - ID.s[i] + 1;
      if (st.minlen && len < st.minlen) continue;
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
      if (st.blocky && !(ID.dz[i] >= 200)) continue;
      if (st.dt && ID.zb[i] == null) continue;
      rows.push(i);
    }
    const f3 = (v) => (v == null ? '' : (v / 1000).toFixed(3));
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
      { key: 'fde', label: 'D+E', num: true, get: (i) => ID.fde[i] / 1000, html: (i) => f3(ID.fde[i]) },
      { key: 'frk', label: 'R+K', num: true, get: (i) => ID.frk[i] / 1000, html: (i) => f3(ID.frk[i]) },
      { key: 'lde', label: 'Longest D/E run', num: true, get: (i) => ID.lde[i] },
      { key: 'nas', label: 'Acidic stretches', num: true, get: (i) => ID.nas[i] },
      { key: 'nbp', label: 'Basic patches', num: true, get: (i) => ID.nbp[i] },
      { key: 'nrg', label: 'RG clusters', num: true, get: (i) => ID.nrg[i] },
      { key: 'mR', label: 'methyl-R', num: true, get: (i) => ID.mR[i] },
      { key: 'ph', label: 'phospho', num: true, get: (i) => ID.ph[i] },
      { key: 'ck2', label: 'CK2', num: true, get: (i) => ID.ck2[i] },
      { key: 'kac', label: 'acyl-K', num: true, get: (i) => ID.kac[i] },
      { key: 'ub', label: 'ubiquitin-K', num: true, get: (i) => ID.ub[i] },
      { key: 'sumo', label: 'SUMO-K', num: true, get: (i) => ID.sumo[i] },
      { key: 'glu', label: 'Glu', num: true, get: (i) => ID.glu[i] },
      { key: 'zb', label: 'dtIDR basic z', num: true, get: (i) => z(ID.zb[i]), html: (i) => zChip(ID.zb[i], 'basic') },
      { key: 'za', label: 'dtIDR acidic z', num: true, get: (i) => z(ID.za[i]), html: (i) => zChip(ID.za[i], 'acidic') },
    ];
    const save = (over) => {
      Object.assign(st, over);
      const p = new URLSearchParams();
      if (st.q) p.set('q', st.q);
      if (st.sort !== 'len') p.set('sort', st.sort);
      if (!st.desc) p.set('asc', '1');
      if (st.page > 1) p.set('page', st.page);
      if (st.minlen) p.set('minlen', st.minlen);
      if (st.term) p.set('term', st.term);
      if (st.ncpr) p.set('ncpr', st.ncpr);
      for (const k of ['as', 'rg', 'bp', 'mixed', 'blocky', 'mr', 'ck2', 'glu', 'dt']) if (st[k]) p.set(k, '1');
      setParams('idrs', p);
    };
    const draw = () => table($('#it'), {
      caption: 'IDRs', columns, rows, sortKey: st.sort, desc: st.desc, page: st.page, csvName: 'idrs.tsv',
      onState: (o) => { save(o); draw(); $('#it').scrollIntoView({ block: 'nearest' }); },
    });
    draw();
    const xs = ID.fcr.map((v) => v / 1000), ys = ID.ncpr.map((v) => v / 1000);
    $('#iplot').innerHTML = '<h3 style="margin:0 0 6px;font-size:1.05rem">Charge of every IDR</h3>' +
      density({ xs, ys, sel: rows.length < n ? rows : null, xLabel: 'fraction of charged residues', yLabel: 'net charge per residue', xr: [0, 1], yr: [-1, 1], bins: 50 }) +
      '<p class="caveat">K and R +1, D and E −1, His 0, as glu.sequence.charge_vector counts them.</p>';
    $('#if').addEventListener('change', (e) => {
      const t = e.target;
      if (t.type === 'checkbox') save({ [t.name]: t.checked, page: 1 });
      else if (t.name === 'minlen') save({ minlen: +t.value || 0, page: 1 });
      else save({ [t.name]: t.value, page: 1 });
      route();
    });
    $('#if').addEventListener('submit', (e) => e.preventDefault());
  }

  // ---------------------------------------------------------------- protein page
  const EMBOSS = { K: 1, R: 1, H: 0.5, D: -1, E: -1 };
  function chargeProfile(seq, bins) {
    const n = seq.length, w = 2;
    const q = Array.from(seq, (c) => EMBOSS[c] || 0);
    const m = new Float32Array(n);
    for (let i = 0; i < n; i++) { let s = 0, k = 0; for (let j = Math.max(0, i - w); j <= Math.min(n - 1, i + w); j++) { s += q[j]; k++; } m[i] = s / k; }
    const out = [];
    for (let b = 0; b < bins; b++) {
      const a = Math.floor(b * n / bins), e = Math.max(a + 1, Math.floor((b + 1) * n / bins));
      let s = 0; for (let i = a; i < e; i++) s += m[i];
      out.push(s / (e - a));
    }
    return out;
  }

  async function proteinView(key, params, token) {
    await loadProteins();
    const K = key.toUpperCase();
    const p = byAcc.has(K) ? byAcc.get(K) : byGene.get(K);
    if (p == null) { view().innerHTML = '<div class="panel">No reviewed human protein called <b>' + esc(key) + '</b>.</div>'; return; }
    const acc = PR.acc[p];
    const [sh] = await Promise.all([loadShard(acc), loadIdrs()]);
    if (token !== routeToken) return;
    const R = sh.prot[acc];
    const names = sh.names;
    const seq = R.seq, n = seq.length;
    const hl = (params.get('f') || '').match(/^(\d+)-(\d+)$/);
    const W = 1000, L = 150, PW = W - L - 10;
    const x = (pos) => L + (pos - 1) / Math.max(1, n - 1) * PW;
    const wd = (s, e) => Math.max(1.5, (e - s + 1) / n * PW);
    let y = 8;
    const rows = [];
    const row = (label, h, inner) => { rows.push('<text class="lab" x="' + (L - 8) + '" y="' + (y + h / 2 + 4) + '" text-anchor="end">' + esc(label) + '</text>' + inner(y, h)); y += h + 8; };
    const seg = (s, e, yy, h, style, title, attrs = '') => '<rect x="' + x(s).toFixed(1) + '" y="' + yy + '" width="' + wd(s, e).toFixed(1) + '" height="' + h + '" rx="2" style="' + style + '" data-s="' + s + '" data-e="' + e + '"' + attrs + '><title>' + esc(title) + '</title></rect>';
    const name = (k) => (names[k] ? names[k][1] : '');
    row('Domains', 16, (yy, h) => (R.dom || []).map((d) => seg(d[0], d[1], yy, h, 'fill:var(--domain);stroke:var(--domain-ink);stroke-width:0.6', name(d[2]) + ' ' + d[0] + '–' + d[1])).join(''));
    if ((R.fam || []).length) row('Families', 8, (yy, h) => R.fam.map((d) => seg(d[0], d[1], yy, h, 'fill:var(--family)', name(d[2]) + ' ' + d[0] + '–' + d[1])).join(''));
    row('Disorder (MobiDB)', 8, (yy, h) => '<line x1="' + L + '" x2="' + (L + PW) + '" y1="' + (yy + h / 2) + '" y2="' + (yy + h / 2) + '" style="stroke:var(--line-strong)"/>' +
      (R.dis || []).map((d) => seg(d[0], d[1], yy, h, 'fill:var(--dis)', 'Disordered ' + d[0] + '–' + d[1])).join(''));
    const flankRow = (label, idx, pctIdx, pole) => row(label, 14, (yy, h) => (R.fl || []).map((f) => {
      const zv = f[idx], pc = f[pctIdx];
      const a = pc == null ? 0.05 : Math.max(0.05, (pc / 1000 - 0.5) * 2);
      return seg(f[1], f[2], yy, h, 'fill:var(--' + pole + ');fill-opacity:' + a.toFixed(2) + ';stroke:var(--' + pole + ');stroke-opacity:0.5;stroke-width:0.5',
        (f[0] ? 'C' : 'N') + '-side flank ' + f[1] + '–' + f[2] + ' (domain ' + f[3] + '–' + f[4] + '): ' + pole + ' z ' + fz(zv) + ', percentile ' + (pc / 10).toFixed(1) + (f[7] >= 0 ? ', ' + CHEM_SHORT[f[7]] : ''));
    }).join(''));
    flankRow('dtIDR, basic pole', 5, 8, 'basic');
    flankRow('dtIDR, acidic pole', 6, 9, 'acidic');
    const bins = Math.min(n, 400), prof = chargeProfile(seq, bins);
    row('Charge (5-aa mean)', 40, (yy, h) => {
      const mid = yy + h / 2, bw = PW / bins;
      let s = '<line x1="' + L + '" x2="' + (L + PW) + '" y1="' + mid + '" y2="' + mid + '" style="stroke:var(--line-strong)"/>';
      prof.forEach((v, b) => {
        if (!v) return;
        const hh = Math.abs(v) * h / 2;
        s += '<rect x="' + (L + b * bw).toFixed(1) + '" y="' + (v > 0 ? mid - hh : mid).toFixed(1) + '" width="' + (bw + 0.2).toFixed(2) + '" height="' + hh.toFixed(1) + '" style="fill:var(--' + (v > 0 ? 'basic' : 'acidic') + ')"/>';
      });
      return s + '<title>EMBOSS charge (K, R +1; H +0.5; D, E −1), sliding 5-residue mean, averaged over ' + (n / bins).toFixed(1) + ' residues per bar</title>';
    });
    row('Acidic stretches', 8, (yy, h) => (R.as || []).map((d) => seg(d[0], d[1], yy, h, 'fill:var(--acidic)', 'Acidic stretch ' + d[0] + '–' + d[1] + ' (glu.sequence.acidic_stretches)')).join(''));
    row('Basic patches', 8, (yy, h) => (R.bp || []).map((d) => seg(d[0], d[1], yy, h, 'fill:var(--basic)', 'Basic patch ' + d[0] + '–' + d[1] + ' (net ≥ +3 in 10 aa)')).join(''));
    row('RG/RGG clusters', 8, (yy, h) => (R.rg || []).map((d) => seg(d[0], d[1], yy, h, 'fill:var(--rg);stroke:var(--rg-ink);stroke-width:0.6', 'RG/RGG cluster ' + d[0] + '–' + d[1])).join(''));
    const ticks = (label, list, color, title, tall = () => false, hollow = null) => row(label, 14, (yy, h) => (list || []).map((v) => {
      const pos = Array.isArray(v) ? v[0] : v, xx = x(pos).toFixed(1), t = tall(v);
      return '<line x1="' + xx + '" x2="' + xx + '" y1="' + (t ? yy - 2 : yy + 3) + '" y2="' + (yy + h) + '" style="stroke:' + color(v) + ';stroke-width:' + (t ? 2.2 : 1.4) + '" data-s="' + pos + '" data-e="' + pos + '"><title>' + esc(title(v)) + '</title></line>' +
        (hollow ? '<circle cx="' + xx + '" cy="' + (yy - 1) + '" r="3" style="fill:' + (hollow(v) ? 'var(--surface)' : color(v)) + ';stroke:' + color(v) + '"/>' : '');
    }).join(''));
    const SRC = (bits) => ['iPTMnet', 'UniProt', 'PTMeXchange'].filter((_, k) => bits & (1 << k)).join(', ');
    const CONF = ['', ', Bronze', ', Silver', ', Gold'];
    const lysTitle = (what) => (v) => 'K' + v[0] + ' ' + what + ' (' + SRC(v[1]) + CONF[v[2] || 0] + ')';
    ticks('methyl-R', R.mR, () => 'var(--meth)', (v) => 'R' + v[0] + (v[1] ? ' ' + v[1] : '') + ' (PTM_Pipeline ' + B.ptm_release + ')');
    ticks('Phospho (CK2 tall)', R.ph, (v) => (v[1] ? 'var(--ck2)' : 'var(--phospho)'), (v) => seq[v[0] - 1] + v[0] + ' phospho' + (v[1] ? ', by CK2' : '') + ' (' + SRC(v[2]) + ')', (v) => !!v[1]);
    ticks('Acyl-K', R.kac, () => 'var(--kac)', lysTitle('acetyl/acyl'), (v) => v[2] === 3);
    if ((R.kme || []).length) ticks('methyl-K', R.kme, () => 'var(--kme)', lysTitle('methyl'));
    if ((R.ub || []).length) ticks('Ubiquitin-K', R.ub, () => 'var(--phospho)', lysTitle('ubiquitin'), (v) => v[2] === 3);
    if ((R.sumo || []).length) ticks('SUMO-K', R.sumo, () => 'var(--kme)', lysTitle('SUMO'), (v) => v[2] === 3);
    if ((R.adpr || []).length) ticks('ADP-ribosyl', R.adpr, () => 'var(--ck2)', (v) => v[1] + v[0] + ' ADP-ribosylation (UniProt)');
    if ((R.cit || []).length) ticks('Citrulline', R.cit, () => 'var(--meth)', (v) => 'R' + v[0] + ' citrullination (UniProt)');
    ticks('Glutamylation', R.glu, () => 'var(--glu)', (v) => 'E' + v[0] + ': ' + B.glu_classes[v[1]] + (v[2] ? ', TTLL4-dependent' : ''), (v) => v[1] === 0, (v) => v[1] > 1 && v[1] < 5);
    const axisY = y;
    let axis = '<line x1="' + L + '" x2="' + (L + PW) + '" y1="' + axisY + '" y2="' + axisY + '" style="stroke:var(--muted)"/>';
    const step = [10, 20, 50, 100, 200, 250, 500, 1000, 2000, 5000].find((s) => n / s <= 10) || 10000;
    for (let t = 0; t <= n; t += step) { const xx = x(Math.max(1, t)); axis += '<line x1="' + xx + '" x2="' + xx + '" y1="' + axisY + '" y2="' + (axisY + 4) + '" style="stroke:var(--muted)"/><text x="' + xx + '" y="' + (axisY + 16) + '" text-anchor="middle" font-size="11">' + Math.max(1, t) + '</text>'; }
    const H = axisY + 24;
    const band = hl ? '<rect x="' + x(+hl[1]).toFixed(1) + '" y="0" width="' + wd(+hl[1], +hl[2]).toFixed(1) + '" height="' + axisY + '" style="fill:var(--accent);fill-opacity:0.12;stroke:var(--accent);stroke-opacity:0.6"/>' : '';
    const svg = '<svg class="map" id="pmap" viewBox="0 0 ' + W + ' ' + H + '" role="img" aria-label="Protein map of ' + esc(geneOf(p)) + '">' + band + rows.join('') + axis + '</svg>';

    const myIdrs = [];
    for (let i = 0; i < ID.p.length; i++) if (ID.p[i] === p) myIdrs.push(i);
    const legend = [['var(--domain)', 'InterPro domain'], ['var(--dis)', 'disordered'], ['var(--basic)', 'basic: dtIDR pole, charge, patch'], ['var(--acidic)', 'acidic: dtIDR pole, charge, stretch'],
      ['var(--rg)', 'RG/RGG cluster'], ['var(--meth)', 'methyl-R'], ['var(--phospho)', 'phospho'], ['var(--ck2)', 'CK2 phospho'], ['var(--kac)', 'acyl-K'], ['var(--glu)', 'glutamylation (filled: localized and TTLL4-dependent, or published)']]
      .concat([['var(--kme)', 'SUMO-K / methyl-K']])
      .map(([c, l]) => '<span><i style="background:' + c + '"></i>' + l + '</span>').join('');
    view().innerHTML =
      '<div class="prot-head"><h2>' + esc(geneOf(p)) + ' <span class="muted" style="font-size:1rem">' + esc(acc) + '</span></h2><p style="margin:0 0 6px">' + esc(PR.name[p]) + ' · ' + num(n) + ' aa · ' + PR.dis[p] + '% disordered</p>' +
      '<div class="links"><a href="' + PTM_SITE + '#' + acc + '">Methylarginine Site Explorer</a><a href="https://www.uniprot.org/uniprotkb/' + acc + '/entry">UniProt</a>' +
      '<a href="https://alphafold.ebi.ac.uk/entry/' + acc + '">AlphaFold</a><a href="https://mobidb.org/' + acc + '">MobiDB</a></div></div>' +
      '<div class="panel"><div class="map-scroll">' + svg + '</div><div class="legend">' + legend + '</div><p class="caveat">Hover a mark for details; click a segment to read its sequence.</p><div id="segment"></div></div>' +
      '<div id="pfl"></div><div id="pidr" style="margin-top:16px"></div>';
    const showSeg = (s, e, label) => {
      s = Math.max(1, s - (s === e ? 15 : 0)); e = Math.min(n, e + (s === e ? 15 : 0));
      const marks = new Set();
      for (const k of ['mR', 'ph', 'kac', 'kme', 'ub', 'sumo', 'glu', 'adpr', 'cit']) for (const v of R[k] || []) marks.add(Array.isArray(v) ? v[0] : v);
      let html = '';
      for (let i = s; i <= e; i++) {
        const c = seq[i - 1];
        const cls = ('KR'.includes(c) ? 'b' : 'DE'.includes(c) ? 'a' : '') + (marks.has(i) ? ' s' : '');
        html += cls ? '<span class="' + cls.trim() + '" title="' + c + i + '">' + c + '</span>' : c;
        if ((i - s + 1) % 10 === 0) html += ' ';
      }
      $('#segment').innerHTML = '<h3 style="margin:10px 0 4px;font-size:1rem">' + esc(label || '') + ' ' + s + '–' + e + '</h3><div class="seq">' + html + '</div>' +
        '<p class="caveat">Basic residues in blue, acidic in red; modified residues underlined.</p>';
    };
    $('#pmap').addEventListener('click', (e) => {
      const t = e.target.closest('[data-s]');
      if (t) showSeg(+t.dataset.s, +t.dataset.e, t.querySelector('title') ? t.querySelector('title').textContent.split(':')[0].split(' (')[0] : '');
    });
    if (hl) showSeg(+hl[1], +hl[2], 'Selected');
    const fl = (R.fl || []).map((f, k) => k);
    table($('#pfl'), {
      caption: 'dtIDR flanks of ' + geneOf(p), rows: fl, sortKey: 'fs', desc: false, page: 1, rowHeader: 0,
      columns: [
        { key: 'fs', label: 'Flank', get: (k) => R.fl[k][1] + '-' + R.fl[k][2], sort: (k) => R.fl[k][1], html: (k) => '<a href="#' + acc + '?f=' + R.fl[k][1] + '-' + R.fl[k][2] + '">' + R.fl[k][1] + '–' + R.fl[k][2] + '</a>' },
        { key: 'side', label: 'Side', get: (k) => (R.fl[k][0] ? 'C' : 'N') },
        { key: 'dom', label: 'Domain', get: (k) => R.fl[k][3] + '-' + R.fl[k][4] },
        { key: 'zb', label: 'Basic z', num: true, get: (k) => z(R.fl[k][5]), html: (k) => zChip(R.fl[k][5], 'basic') },
        { key: 'pb', label: 'Basic percentile', num: true, get: (k) => R.fl[k][8] / 10 },
        { key: 'za', label: 'Acidic z', num: true, get: (k) => z(R.fl[k][6]), html: (k) => zChip(R.fl[k][6], 'acidic') },
        { key: 'pa', label: 'Acidic percentile', num: true, get: (k) => R.fl[k][9] / 10 },
        { key: 'chem', label: 'Chemistry', get: (k) => CHEM_SHORT[R.fl[k][7]] || '' },
      ],
    });
    const f3 = (v) => (v / 1000).toFixed(3);
    table($('#pidr'), {
      caption: 'IDRs of ' + geneOf(p) + ' (≥20 aa)', rows: myIdrs, sortKey: 's', desc: false, page: 1,
      columns: [
        { key: 's', label: 'IDR', get: (i) => ID.s[i] + '-' + ID.e[i], sort: (i) => ID.s[i], html: (i) => '<a href="#' + acc + '?f=' + ID.s[i] + '-' + ID.e[i] + '">' + ID.s[i] + '–' + ID.e[i] + '</a>' },
        { key: 'len', label: 'Length', num: true, get: (i) => ID.e[i] - ID.s[i] + 1 },
        { key: 'ncpr', label: 'NCPR', num: true, get: (i) => f3(ID.ncpr[i]) },
        { key: 'fcr', label: 'FCR', num: true, get: (i) => f3(ID.fcr[i]) },
        { key: 'dz', label: 'δ z (blocks)', num: true, get: (i) => z(ID.dz[i]), html: (i) => fz(ID.dz[i]) },
        { key: 'lde', label: 'Longest D/E run', num: true, get: (i) => ID.lde[i] },
        { key: 'nas', label: 'Acidic stretches', num: true, get: (i) => ID.nas[i] },
        { key: 'nbp', label: 'Basic patches', num: true, get: (i) => ID.nbp[i] },
        { key: 'nrg', label: 'RG clusters', num: true, get: (i) => ID.nrg[i] },
        { key: 'mR', label: 'methyl-R', num: true, get: (i) => ID.mR[i] },
        { key: 'ph', label: 'phospho', num: true, get: (i) => ID.ph[i] },
        { key: 'glu', label: 'Glu', num: true, get: (i) => ID.glu[i] },
      ],
    });
    document.title = geneOf(p) + ' · IDR Explorer';
  }

  // ---------------------------------------------------------------- sources
  async function aboutView() {
    const c = B.counts;
    const inputs = Object.entries(B.inputs).map(([k, v]) => '<tr><td>' + esc(k) + '</td><td>' + esc(v.repo) + '</td><td><code>' + esc(v.path) + '</code></td><td><code>' + esc(v.md5.slice(0, 10)) + '…</code></td></tr>').join('');
    const commit = (repo, sha) => sha ? '<a href="https://github.com/Shechterlab/' + repo + '/commit/' + sha + '"><code>' + sha.slice(0, 9) + '</code></a>' : 'not recorded';
    view().innerHTML = '<div class="about"><h2>Sources and rules</h2>' +
      '<div class="panel"><p>This site is generated, not edited: <code>python3 Analysis/scripts/build_idr_explorer.py</code> in a MicroDomains clone with PTM_Pipeline, glu_proteomics and this repository cloned beside it. It reuses each repository’s own code for what that repository owns, and its sequence definitions are <a href="' + GH.site + 'idrkit">idrkit</a>, the package in this repository that glu_proteomics imports too, so a definition changes in one place.</p>' +
      '<dl><dt>Built</dt><dd>' + esc(B.generated) + ' from PTM_Pipeline ' + commit('PTM_Pipeline', B.commits.PTM_Pipeline) + ' (release ' + esc(B.ptm_release) + '), idr-explorer ' + commit('idr-explorer', B.commits['idr-explorer']) + ', MicroDomains <code>' + esc((B.commits.MicroDomains || '').slice(0, 9)) + '</code> and glu_proteomics <code>' + esc((B.commits.glu_proteomics || '').slice(0, 9)) + '</code> (private).</dd>' +
      '<dt>dtIDR scores</dt><dd>MicroDomains: basic score <code>dtidr_score_z</code> and the acidic score on the same 126,706 windows (acidic_dtidr_score.py), chemistry classes. Shechter Lab, unpublished.</dd>' +
      '<dt>IDR</dt><dd>' + esc(B.rules.idr) + '.</dd>' +
      '<dt>Acidic stretch</dt><dd>glu_proteomics <code>glu.sequence.acidic_stretches</code>: a ' + B.rules.acidic_stretch.window + '-residue window with ≥ ' + B.rules.acidic_stretch.min_acidic + ' D/E and ≤ ' + B.rules.acidic_stretch.max_basic + ' K/R, merged across gaps ≤ ' + B.rules.acidic_stretch.merge_gap + ', trimmed to D/E, ≥ ' + B.rules.acidic_stretch.min_length + ' residues.</dd>' +
      '<dt>RG/RGG cluster</dt><dd><code>glu.sequence.rg_clusters</code>: ≥ ' + B.rules.rg_cluster.min_rg + ' RG dipeptides within ' + B.rules.rg_cluster.window + ' residues.</dd>' +
      '<dt>Basic patch</dt><dd>' + esc(B.rules.basic_patch) + '.</dd>' +
      '<dt>Domains</dt><dd>InterPro 108.0 rows exactly as the Methylarginine Site Explorer draws them (PTM_Pipeline <code>scripts/36c</code>).</dd>' +
      '<dt>Methylarginines</dt><dd>PTM_Pipeline release ' + esc(B.ptm_release) + ' (' + num(c.methylarginines) + ' sites): DeAngelo, Silverstein &amp; Shechter, Cell Chem Biol 2026, doi:10.1016/j.chembiol.2026.09.011. CC BY-NC-SA 4.0; cite the studies in the release’s studies table.</dd>' +
      '<dt>Phosphorylation</dt><dd>' + num(c.phosphosites) + ' sites (' + num(c.ck2_phosphosites) + ' attributed to CK2) from PTM_Pipeline’s iPTMnet table (Huang et al., Nucleic Acids Res 46:D542, 2018; CC BY-NC-SA 4.0; no PhosphoSitePlus) and UniProt 2026_03 MOD_RES (CC BY 4.0), the table the analyses use (MicroDomains export_ptm_sites.py).</dd>' +
      '<dt>Lysine marks</dt><dd>Acetyl/acyl (' + num(c.acyl_lysines) + '), ubiquitin (' + num(c.ubiquitin_lysines) + '), SUMO (' + num(c.sumo_lysines) + ') and methyl (' + num(c.methyl_lysines) + ') lysines from iPTMnet, UniProt and PTMeXchange’s reprocessed human builds (PRIDE PXD079218, PXD068989, PXD068489, PXD076372; CC0), each PTMeXchange site graded Gold, Silver or Bronze and placed on the 2026_03 sequence by its own peptides (PTM_Pipeline scripts/54).</dd>' +
      '<dt>Other marks</dt><dd>ADP-ribosylation (' + num(c.adp_ribosyl) + ') and citrullination (' + num(c.citrulline) + ') from UniProt 2026_03 (PTM_Pipeline scripts/53; CC BY 4.0).</dd>' +
      '<dt>Glutamylation</dt><dd>The TTLL4 site list of a re-analysis of PRIDE PXD065251 (Dey et al., J Cell Sci 2026, doi:10.1242/jcs.264734), its known sites (NPM1 E126) and UniProt’s polyglutamylation annotations: ' + num(c.glu_sites) + ' sites, ' + num(c.glu_sites_dropped) + ' not placed on a glutamate. Classes: ' + B.glu_classes.map(esc).join('; ') + '.</dd>' +
      '<dt>Patterning</dt><dd>' + esc(B.rules.patterning) + '. δ matches localCIDER’s (Das &amp; Pappu, PNAS 110:13392, 2013) to 1e-16.</dd>' +
      '<dt>Reference, disorder</dt><dd>UniProt 2026_03 reviewed human (CC BY 4.0); MobiDB (Piovesan et al., Nucleic Acids Res 53:D495, 2025; CC BY 4.0 as listed); InterPro 108.0 (CC0).</dd>' +
      '<dt>Analyses</dt><dd>The dtIDR scores, the statistics on the Ideas page and their figures are the Shechter lab’s MicroDomains analyses (unpublished; the repository is private until publication), built from the public inputs above.</dd>' +
      '<dt>Not in the reference</dt><dd>' + (c.flank_proteins_not_in_reference || []).map(esc).join(', ') + ' (flank accessions without a 2026_03 sequence, left out).</dd></dl></div>' +
      '<div class="panel"><h3 style="margin-top:0">Inputs</h3><div class="table-wrapper" tabindex="0" role="region" aria-label="Inputs and checksums"><table class="data"><thead><tr><th>Input</th><th>Repository</th><th>Path</th><th>MD5</th></tr></thead><tbody>' + inputs + '</tbody></table></div></div>' +
      '<div class="panel"><h3 style="margin-top:0">Read before quoting</h3><ul>' +
      '<li>Modification sites are reported sites: they track study effort as well as biology.</li>' +
      '<li>Scores live on a lattice (most windows are 50 aa; disorder fractions are often exactly 0 or 1); ties are expected.</li>' +
      '<li>The acidic pole is scored on the full proteome since 27 October 2026; two earlier claims (the histone fold depleted at the acidic pole, histone-like proteins depleted) did not survive the move.</li>' +
      '<li>Glutamylation sites from MS favour peptides beside K/R; most are region-level.</li></ul></div></div>';
  }

  // ---------------------------------------------------------------- search and theme
  const input = $('#q'), sug = $('#suggest');
  const matches = (t) => {
    const T = t.trim().toUpperCase();
    if (!T || !PR) return [];
    const out = [];
    if (byAcc.has(T)) out.push(byAcc.get(T));
    // genes starting with the text, the exact and shortest first (NAP1L1 before NAP1L6P), then accessions
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
  $('#theme-btn').addEventListener('click', () => {
    const cur = document.documentElement.getAttribute('data-theme') || (matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light');
    const next = cur === 'dark' ? 'light' : 'dark';
    document.documentElement.setAttribute('data-theme', next);
    try { localStorage.setItem('idr-theme', next); } catch (e) { /* private mode */ }
  });

  loadBuild().then(() => { $('#built').textContent = 'Data built ' + B.generated + '.'; }).catch(() => {});
  route();
})();
