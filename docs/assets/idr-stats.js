/* IDR Explorer: the Overview's tests, for the whole dataset or a selection.

   The same tests as Analysis/scripts/idr_overview_stats.py (Shechterlab/MicroDomains), on the same per-protein counts
   (data/residues.json, written from idr_universe.background) and sites (data/sites.json), run with the Methylarginine
   Site Explorer's own implementations (ptm-core.js: mhTest, fisherP, hyperSf, bhAdjust, mwuP, erfc), so a selection is
   tested exactly as the whole dataset was; tests/overview.test.mjs checks the whole-dataset results against the Python
   tables. Pure functions over the data files, no DOM.

   A selection is { prots: Set of protein indices or null (every protein), keep: (site index) -> bool or null (every
   site) }. Sites a selection does not keep are left out of both sides of a residue test, as the Methylarginine Site
   Explorer leaves out the release's unselected sites. */
window.IDRStats = (function () {
  'use strict';

  function make(C, D) {
    const OV = D.overview, S = D.sites, RS = D.residues, P = D.proteins, I = D.idrs;
    const CTX = OV.contexts.map((c) => c[0]);
    const KINDS = OV.kinds.map((k) => k[0]);
    const KRES = Object.fromEntries(OV.kinds.map((k) => [k[0], k[2]]));
    const GI = Object.fromEntries(RS.groups.map((g, i) => [g, i]));
    const n = S.p.length;
    // a test needs enough information (idr_universe.MIN_IN_CONTEXT residues in the context, MIN_STRATA informative proteins)
    const MIN_IN = OV.settings.min_in_context, MIN_ST = OV.settings.min_strata;
    // the sites of each kind (tested ones: a residue the kind can carry)
    const ofKind = KINDS.map(() => []);
    for (let i = 0; i < n; i++) if (S.t[i]) ofKind[S.k[i]].push(i);
    const keepAll = (sel) => !sel || (!sel.prots && !sel.keep);
    const inSel = (sel, i) => !sel || ((!sel.prots || sel.prots.has(S.p[i])) && (!sel.keep || sel.keep(i)));

    // Sites of `kind` in context `ctx` against the unmodified residues of the same type in the same proteins: pooled
    // (Haldane-Anscombe odds ratio, Fisher's exact test) and within protein (Mantel-Haenszel, each protein a stratum).
    // idr_universe.residue_test.
    function residueTest(kind, ctx, sel) {
      const k = KINDS.indexOf(kind), ci = CTX.indexOf(ctx), gi = GI[KRES[kind]];
      const per = new Map();   // protein -> [selected in ctx, selected, known in ctx, known]
      for (const i of ofKind[k]) {
        const p = S.p[i];
        let r = per.get(p);
        if (!r) per.set(p, (r = [0, 0, 0, 0]));
        const inn = (S.ctx[i] >> ci) & 1;
        r[2] += inn; r[3]++;
        if (inSel(sel, i)) { r[0] += inn; r[1]++; }
      }
      const tables = [];
      let a = 0, b = 0, c = 0, d = 0;
      for (const [p, r] of per) {
        if (!r[1]) continue;
        const rec = RS.bg[p][gi];
        if (!rec) continue;
        const tot = rec[0], cin = rec[1 + ci];
        const t = [r[0], r[1] - r[0], cin - r[2], (tot - cin) - (r[3] - r[2])];
        tables.push(t);
        a += t[0]; b += t[1]; c += t[2]; d += t[3];
      }
      const m = tables.length ? C.mhTest(tables) : { or: null, lo: null, hi: null, p: null, used: 0 };
      return { kind, context: ctx, sites_in: a, sites: a + b, other_in: c, other: c + d, proteins: tables.length, tested: a + c >= MIN_IN && m.used >= MIN_ST && m.p != null,
        pooled_or: tables.length ? ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5)) : null, pooled_p: tables.length ? C.fisherP(a, b, c, d) : null,
        mh_or: m.or, mh_lo: m.lo, mh_hi: m.hi, mh_p: m.p, mh_strata: m.used };
    }

    // Observed phospho-S/T in context `ctx` against the other phospho-S/T of the same proteins: the share that are motif
    // matches of kinase row `row` (bit of sites.km). idr_overview_stats.kinase_groups.
    const ph = KINDS.indexOf('ph');
    function kinaseTest(ctx, row, sel) {
      const ci = CTX.indexOf(ctx), bit = 1 << row;
      const per = new Map();
      for (const i of ofKind[ph]) {
        if (S.km[i] < 0 || !inSel(sel, i)) continue;
        const p = S.p[i];
        let t = per.get(p);
        if (!t) per.set(p, (t = [0, 0, 0, 0]));
        const inn = (S.ctx[i] >> ci) & 1, hit = S.km[i] & bit ? 1 : 0;
        t[inn ? (hit ? 0 : 1) : (hit ? 2 : 3)]++;
      }
      const tables = [...per.values()];
      const m = C.mhTest(tables);
      const sum = (j) => tables.reduce((s, t) => s + t[j], 0);
      return { row, context: ctx, in_context: sum(0) + sum(1), in_context_match: sum(0), outside: sum(2) + sum(3), outside_match: sum(2), tested: sum(0) + sum(1) >= MIN_IN && m.used >= MIN_ST && m.p != null,
        mh_or: m.or, mh_lo: m.lo, mh_hi: m.hi, mh_p: m.p, mh_strata: m.used };
    }

    // Protein classes among a set of proteins (indices) against all reviewed proteins: one-sided hypergeometric, fold
    // enrichment, BH across the classes. idr_overview_stats.protein_sets.
    function proteinSets(members) {
      const N = P.acc.length, n = members.length;
      const rows = OV.classes.map(([key], b) => {
        let M = 0, k = 0;
        for (let i = 0; i < N; i++) if (P.cls[i] & (1 << b)) M++;
        for (const i of members) if (P.cls[i] & (1 << b)) k++;
        return { class: key, k, n, M, N, fold: M && n ? (k / n) / (M / N) : null, p: k ? C.hyperSf(k, N, M, n) : 1 };
      });
      C.bhAdjust(rows, (r) => r.p);
      return rows;
    }
    const hasProteins = (bit) => { const out = []; for (let i = 0; i < P.acc.length; i++) if (P.has[i] & (1 << bit)) out.push(i); return out; };

    // Occupancy of sites (sites.occ: phospho, % of the protein; sites.ac: acetyl, fraction) inside a context against
    // outside it: quantiles and Mann-Whitney. idr_overview_stats.stoichiometry.
    const quant = (v) => {
      if (!v.length) return null;
      const s = v.slice().sort((x, y) => x - y), q = (f) => { const h = (s.length - 1) * f, lo = Math.floor(h); return s[lo] + (s[Math.min(lo + 1, s.length - 1)] - s[lo]) * (h - lo); };
      return [q(0.05), q(0.25), q(0.5), q(0.75), q(0.95), s.length];
    };
    function stoich(which, ctx, sel) {
      const ci = CTX.indexOf(ctx), pairs = which === 'occ' ? S.occ : S.ac, scale = which === 'occ' ? 0.01 : 1;
      const x = [], y = [];
      for (const [i, v] of pairs) {
        if (!inSel(sel, i)) continue;
        ((S.ctx[i] >> ci) & 1 ? x : y).push(v * scale);
      }
      const [p, mx, my] = C.mwuP(x, y);
      return { context: ctx, n_in: x.length, n_out: y.length, q_in: quant(x), q_out: quant(y), median_in: mx, median_out: my, mwu_p: p };
    }

    // The charge change observed modifications make in each IDR, against random placement of as many residues of the
    // same type in the same IDR (idrs.json dd_<kind>, rd_<kind>): a sign test on the IDRs where the two differ.
    // idr_overview_stats.charge_shift_summary.
    function signP(above, below) {
      const m = above + below;
      if (!m) return null;
      const z = Math.max(0, Math.abs(above - m / 2) - 0.5) / Math.sqrt(m / 4);
      return Math.min(1, C.erfc(z / Math.SQRT2));
    }
    function chargeShift(kind, prots) {
      const dd = I['dd_' + kind], rd = I['rd_' + kind];
      let idrs = 0, above = 0, below = 0, up = 0;
      for (let j = 0; j < dd.length; j++) {
        if (dd[j] == null || rd[j] == null || (prots && !prots.has(I.p[j]))) continue;
        idrs++;
        if (dd[j] > 0) up++;
        if (dd[j] > rd[j]) above++; else if (dd[j] < rd[j]) below++;
      }
      return { kind, idrs, above, below, sign_p: signP(above, below), share_delta_up: idrs ? up / idrs : null };
    }

    return { residueTest, kinaseTest, proteinSets, hasProteins, stoich, chargeShift, signP, keepAll, inSel, CTX, KINDS, KRES, ofKind };
  }
  return { make };
})();
