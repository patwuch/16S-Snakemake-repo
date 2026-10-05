#!/usr/bin/env python3
"""Check the factorial report's Beta diversity section against its tables.

check(html) compares, cell by cell, the PERMANOVA table (q, floor tag, verdict, R²)
and the PERMDISP table (q, verdict) with tables/diversity_factorial_Day17.tsv, the
output of run_diversity_factorial.R, and recomputes every number the beta prose,
captions, findings and summary quote: the PCoA axis shares (recomputed on the
Day17 mice, as plot_diversity_factorial.py draws them), the R² ranges and the
residual, and the PERMDISP q and raw p values. Returns a list of misses.

A displayed number passes when it equals the table value rounded to the number's
own displayed precision, so the check does not depend on a formatting rule.

Called by check_report_factorial.py, inside the qiime2 container (the PCoA needs
qiime tools export and scikit-bio).
"""
import os
import re
import warnings

import numpy as np
import pandas as pd

import plot_diversity_factorial as pdf

TERMS = {"Treatment": "Treatment", "Radiation": "Radiation", "Interaction": "Treatment:Radiation"}
METRICS = {"Jaccard": ("jaccard", "any"), "Bray-Curtis": ("bray_curtis", "any"),
           "Unweighted UniFrac · Greengenes2": ("unweighted_unifrac", "Greengenes2"),
           "Weighted UniFrac · Greengenes2": ("weighted_unifrac", "Greengenes2")}
FLOOR_P = 1 / 10_000  # (0 hits + 1) / (9,999 + 1)


def shown_equals(shown, value):
    """'0.0002' equals 0.00018…; '0.70' equals 0.6959; '16%' equals 0.163 (as a %)."""
    s = shown.strip().rstrip("%")
    decimals = len(s.split(".")[1]) if "." in s else 0
    return abs(float(s) - value) <= 0.5 * 10 ** -decimals + 1e-12


def table_rows(section, header):
    """(group, metric, cells…) for each row of the table whose header row contains `header`."""
    i = section.index(header)
    body = section[section.index("<tbody>", i):section.index("</tbody>", i)]
    rows, group = [], None
    for r in re.findall(r"<tr[^>]*>(.*?)</tr>", body, re.S):
        c = [re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", x)).strip()
             for x in re.findall(r"<td[^>]*>(.*?)</td>", r, re.S)]
        raw = re.findall(r"<td[^>]*>(.*?)</td>", r, re.S)
        if 'rowspan' in r.split("</td>")[0]:
            group, c, raw = c[0], c[1:], raw[1:]
        rows.append((group, c, raw))
    return rows


_SHARES = {}


def pcoa_shares(meta):
    """Share of variation on the first two PCoA axes per metric, Day17 mice only.
    Cached per process: each metric needs a qiime export, which is slow."""
    key = tuple(meta.index)
    if key in _SHARES:
        return _SHARES[key]
    from skbio import DistanceMatrix
    from skbio.stats.ordination import pcoa
    out = {}
    warnings.filterwarnings("ignore", message="The result contains negative eigenvalues")
    for metric, db, fname, _label in pdf.BETA_PANELS:
        dm = pd.read_csv(pdf.export_qza(os.path.join(pdf.CORE, fname), "distance-matrix.tsv"), sep="\t", index_col=0)
        ids = [s for s in meta.index if s in dm.index]
        out[metric] = float(pcoa(DistanceMatrix(dm.loc[ids, ids].values, ids=ids)).proportion_explained[:2].sum())
    _SHARES[key] = out
    return out


def check(html):
    misses = []
    b0, b1 = html.index("<h2>Beta diversity</h2>"), html.index("<h2>Differential abundance</h2>")
    beta = html[b0:b1]
    text = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", html))
    btext = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", beta))

    def expect(value, why, where=text):
        if value not in where:
            misses.append(f"beta {why}: '{value}' not in report")

    st = pd.read_csv(os.path.join(pdf.TAB_DIR, "diversity_factorial_Day17.tsv"), sep="\t")
    st = st[st.db != "Silva138"]
    get = lambda analysis, metric, term: st[(st.analysis == analysis) & (st.metric == metric[0])
                                           & (st.db == metric[1]) & (st.term == term)].iloc[0]

    # ---- PERMANOVA table, cell by cell ----
    perm = table_rows(beta, "<th>Term</th>")
    if len(perm) != 12:
        misses.append(f"beta PERMANOVA table has {len(perm)} rows, expected 12")
    for term, (metric, q_txt, verdict, r2_txt), raw in ((g, c, r) for g, c, r in perm):
        r = get("beta_permanova", METRICS[metric], TERMS[term])
        q_num = q_txt.replace("(floor)", "").strip()
        if not shown_equals(q_num, r.p_BH):
            misses.append(f"beta PERMANOVA {term}/{metric}: q shown {q_num}, table {r.p_BH:.4g}")
        if ("(floor)" in q_txt) != (r.p_raw <= FLOOR_P):
            misses.append(f"beta PERMANOVA {term}/{metric}: floor tag {'shown' if '(floor)' in q_txt else 'missing'}, p = {r.p_raw}")
        if ("<strong>" in raw[1]) != (r.p_BH < 0.05):
            misses.append(f"beta PERMANOVA {term}/{metric}: q bolding does not match q < 0.05")
        if verdict.lower().startswith("significant") != (r.p_BH < 0.05):
            misses.append(f"beta PERMANOVA {term}/{metric}: verdict '{verdict}' but q = {r.p_BH:.3g}")
        if not shown_equals(r2_txt, r.effect * 100):
            misses.append(f"beta PERMANOVA {term}/{metric}: R² shown {r2_txt}, table {r.effect * 100:.1f}%")

    # ---- PERMDISP table: Cell is betadisper on the six groups; the term rows split it by factor ----
    disp = table_rows(beta, "<th>Grouping</th>")
    if len(disp) != 16:
        misses.append(f"beta PERMDISP table has {len(disp)} rows, expected 16")
    for grouping, (metric, q_txt, verdict), raw in ((g, c, r) for g, c, r in disp):
        if grouping.startswith("Cell"):
            r = get("beta_dispersion", METRICS[metric], "Cell")
        else:
            r = get("beta_dispersion_factorial", METRICS[metric], TERMS[grouping])
        if not shown_equals(q_txt, r.p_BH):
            misses.append(f"beta PERMDISP {grouping}/{metric}: q shown {q_txt}, table {r.p_BH:.4g}")
        want = "spread differs" if r.p_BH < 0.05 else ("leans, not significant" if r.p_BH < 0.1 else "equal spread")
        if verdict != want:
            misses.append(f"beta PERMDISP {grouping}/{metric}: verdict '{verdict}', q = {r.p_BH:.3g} gives '{want}'")
        if ("<strong>" in raw[1]) != (r.p_BH < 0.05):
            misses.append(f"beta PERMDISP {grouping}/{metric}: q bolding does not match q < 0.05")

    # ---- PCoA axis shares ----
    meta = pd.read_csv(os.path.join(pdf.STUDY_DIR, "metadata.tsv"), sep="\t", index_col=0,
                       keep_default_na=False, na_values=[""])
    meta = meta[meta.Timepoint == "Day17"]
    sh = {k: round(v * 100) for k, v in pcoa_shares(meta).items()}
    expect(f"They carry {sh['jaccard']}% of the variation for Jaccard, {sh['bray_curtis']}% for Bray-Curtis, "
           f"{sh['unweighted_unifrac']}% for unweighted UniFrac and {sh['weighted_unifrac']}% for weighted UniFrac",
           "PCoA axis shares")

    # ---- R² ranges and residual (R² figure caption, findings, summary) ----
    pm = st[st.analysis == "beta_permanova"]
    pct = lambda s: s * 100
    rng = lambda s: f"{pct(s.min()):.0f}–{pct(s.max()):.0f}%"
    tr, ra, ia = (pm[pm.term == t].set_index("metric").effect for t in ("Treatment", "Radiation", "Treatment:Radiation"))
    resid = 1 - pm.groupby("metric").effect.sum()
    sig_ia = pm[(pm.term == "Treatment:Radiation") & (pm.p_BH < 0.05)]
    ns_ia = pm[(pm.term == "Treatment:Radiation") & (pm.p_BH >= 0.05)]
    if len(sig_ia) != 3 or list(ns_ia.metric) != ["weighted_unifrac"]:
        misses.append("beta: text says the interaction passes on three metrics and fails only on weighted UniFrac")
    expect(f"Treatment is the largest effect everywhere ({rng(tr)}), radiation is {rng(ra).replace('%', '')}%", "R² caption ranges")
    expect(f"the interaction is {rng(sig_ia.effect)} on the three metrics where it is significant but "
           f"{pct(ns_ia.effect.iloc[0]):.0f}% on weighted UniFrac", "R² caption interaction")
    expect(f"The residual, {rng(resid)}, is variation", "R² caption residual")
    if not (tr > ra.reindex(tr.index)).all() or not (tr >= pm.groupby("metric").effect.max().reindex(tr.index)).all():
        misses.append("beta: text says treatment is the largest effect on every metric")
    expect(f"(R² {rng(tr)} against {rng(ra)})", "treatment vs radiation R²")
    expect(f"(q ≤ {sig_ia.p_BH.max():.3f}, R² {rng(sig_ia.effect)})", "interaction finding q and R²")
    expect(f"(q = {ns_ia.p_BH.iloc[0]:.2f}, R² {pct(ns_ia.effect.iloc[0]):.0f}%)", "weighted UniFrac interaction")
    tr_ra = pm[pm.term.isin(["Treatment", "Radiation"])]
    if not (tr_ra.p_BH < 0.05).all():
        misses.append("beta: text says all eight Treatment and Radiation rows clear the bar")
    n_floor = int((tr_ra.p_raw <= FLOOR_P).sum())
    expect(f"{['no', 'one', 'two', 'three', 'four', 'five'][n_floor]} of them at the permutation floor", "floor count")
    wr = get("beta_permanova", METRICS["Weighted UniFrac · Greengenes2"], "Radiation")
    expect(f"Weighted UniFrac's radiation result (q = {wr.p_BH:.3f})", "weighted UniFrac radiation q")
    expect(f"absent on weighted UniFrac (q = {ns_ia.p_BH.iloc[0]:.2f})", "weighted UniFrac interaction q")
    if not (sig_ia.p_BH < 0.01).all():
        misses.append("summary: text says the presence-sensitive metrics catch the interaction at q < 0.01")

    # ---- PERMDISP prose ----
    cell = st[(st.analysis == "beta_dispersion") & (st.term == "Cell")]
    if list(cell[cell.p_BH < 0.05].metric) != ["jaccard"]:
        misses.append("beta: text says Jaccard is the one metric whose spread differs across the six groups")
    expect(f"spread differs across the six groups (q = {cell[cell.metric == 'jaccard'].p_BH.iloc[0]:.3f})", "Jaccard cell q")
    dfac = st[st.analysis == "beta_dispersion_factorial"]
    jt = dfac[(dfac.metric == "jaccard") & (dfac.term == "Treatment")].iloc[0]
    bt = dfac[(dfac.metric == "bray_curtis") & (dfac.term == "Treatment")].iloc[0]
    expect(f"(q = {jt.p_BH:.3f} for Treatment)", "Jaccard treatment dispersion q")
    expect(f"(raw p = {bt.p_raw:.3f}, q = {bt.p_BH:.3f})", "Bray-Curtis treatment dispersion")
    rad_int = dfac[dfac.term.isin(["Radiation", "Treatment:Radiation"])]
    if len(set(rad_int.p_BH.round(2))) != 1:
        misses.append("beta: text quotes one q for every radiation and interaction spread test")
    expect(f"(q = {rad_int.p_BH.iloc[0]:.2f}), so those PERMANOVA results read as shifts", "radiation/interaction spread q")
    pooled = get("beta_dispersion", METRICS["Unweighted UniFrac · Greengenes2"], "Radiation")
    expect(f"near-miss (raw p = {pooled.p_raw:.3f})", "pooled unweighted UniFrac radiation near-miss")
    # Five-ways row: dispersion-clean metrics are those whose six-group spread is equal.
    clean = set(cell[cell.p_BH >= 0.05].metric)
    n_clean_sig = int(sig_ia.metric.isin(clean).sum())
    expect(f"Significant on {n_clean_sig} of {len(clean)} dispersion-clean metrics", "five-ways PERMANOVA row")
    return misses


if __name__ == "__main__":
    import sys
    m = check(open(os.path.join(pdf.OUT_DIR, "report.html")).read())
    print("\n".join(m) if m else "beta checks pass")
    sys.exit(1 if m else 0)
