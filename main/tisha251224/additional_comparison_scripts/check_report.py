#!/usr/bin/env python3
"""Check comparison-Diet-x-Condition-factorial/report.html against its tables.

Every figure the page references must exist and be non-empty, and each number
below is recomputed from its table and must appear in the report text as
formatted. Exits non-zero listing every miss.

Usage, inside the qiime2 container:
    python3 check_report.py
"""
import os
import re
import sys

import numpy as np
import pandas as pd

from da_rows import rows
from plot_interaction_forest import PANELS, joint_q, load_tables

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(THIS_DIR), "comparison-Diet-x-Condition-factorial")
TAB = os.path.join(OUT, "tables")
LN2 = np.log(2)

html = open(os.path.join(OUT, "report.html")).read()
text = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", html))
misses = []


def expect(value, why):
    if value not in text:
        misses.append(f"{why}: '{value}' not in report")


for src in set(re.findall(r'src="([^"]+)"', html)):
    path = os.path.join(OUT, src)
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        misses.append(f"figure missing or empty: {src}")
n_fig = len(set(re.findall(r'src="(figures/[^"]+\.png)"', html)))
expect(f"{n_fig} figures, PNG and SVG", "figure count")

# ---- design ----
meta = pd.read_csv(os.path.join(OUT, "metadata_factorial.tsv"), sep="\t")
expect(f"Samples {len(meta)}", "sample count")

# ---- diversity ----
div = pd.read_csv(os.path.join(TAB, "diversity_factorial.tsv"), sep="\t")
alpha = div[div.analysis == "alpha"]
for r in alpha.itertuples():
    q = r.p_BH
    expect(f"{q:.1e}" if q < 1e-3 else f"{q:.2g}", f"alpha {r.metric} {r.term} q")

perm = div[div.analysis == "beta_permanova"].set_index(["metric", "term"])
for (metric, term), r in perm.iterrows():
    if term == "Residual":
        continue
    q = r.p_BH
    expect("0.0001" if q <= 1e-4 else (f"{q:.4f}" if abs(q - 0.05) < 0.001 else f"{q:.2g}"),
           f"PERMANOVA {metric} {term} q")
    expect(f"{r.effect * 100:.0f}%", f"PERMANOVA {metric} {term} R2")
lo = lambda t: perm.xs(t, level="term").effect
expect(f"{lo('Diet').min() * 100:.0f}–{lo('Diet').max() * 100:.0f}%", "Diet R2 range")
expect(f"{lo('Residual').min() * 100:.0f}–{lo('Residual').max() * 100:.0f}%", "Residual R2 range")

for a in ("beta_dispersion_factorial",):
    for r in div[div.analysis == a].itertuples():
        q = r.p_BH
        expect(f"{q:.4f}" if q < 0.01 else (f"{q:.3f}" if q < 0.1 else f"{q:.2f}"), f"{a} {r.metric} {r.term} q")
for r in div[(div.analysis == "beta_dispersion") & (div.term == "Cell")].itertuples():
    q = r.p_BH
    expect(f"{q:.4f}" if q < 0.01 else f"{q:.2f}", f"PERMDISP cell {r.metric} q")

# ---- differential abundance ----
main = pd.read_csv(os.path.join(TAB, "ANCOMBC2_factorial_Greengenes2.tsv"), sep="\t")
untestable = pd.read_csv(os.path.join(TAB, "ANCOMBC2_factorial_Greengenes2_untestable_absent_from_a_cell.tsv"), sep="\t")
named = main.taxon.map(lambda t: bool(re.match(r"s__.", t.split(";")[-1].strip())))
expect(f"Of 188 features, {len(main) + len(untestable)} are seen in at least six of the 19 samples", "DA taxon set")
expect(f"{len(main)} of them are seen in at least one sample of every one of the six groups", "tested count")
expect(f"{named.sum()} carry a species name, and the other {(~named).sum()}", "species split")
expect(f"The other {len(untestable)} are missing from at least one whole group", "untestable count")
expect(f"Another {len(untestable)} features", "untestable count (note)")
rsir = untestable.columns[untestable.columns.str.startswith("RS+IR")][0]
rs = untestable.columns[untestable.columns.str.startswith("RS (")][0]
absent = untestable[untestable[rsir] == 0]
expect(f"{len(absent)} of the {len(untestable)} are absent from all three RS+IR", "absent from RS+IR")
expect(f"{int((absent[rs] > 0).sum())} of those {len(absent)} are present in unirradiated RS", "present in RS")

_da0 = html.index("Largest up / down (log2FC)")
_da_html = html[_da0:html.index("</table>", _da0)]
for r in rows():
    cells = re.search(rf"(?:<td>|</span>){re.escape(r['label'])}</td><td class=\"n\">(\d+)</td><td class=\"n\">([^<]+)</td><td class=\"n\">([^<]+)</td>", _da_html)
    if not cells:
        misses.append(f"DA row {r['label']}: not found in table")
        continue
    up = f"+{r['max_up']:.2f}" if r["max_up"] is not None else "—"
    down = f"−{r['max_down']:.2f}" if r["max_down"] is not None else "—"
    want = (str(r["n"]), f"{up} / {down}", f"{r['min_q']:.1e}") if r["n"] else ("0", "—", "—")
    got = (cells[1], cells[2], cells[3])
    if want != got:
        misses.append(f"DA row {r['label']}: table has {got}, tables give {want}")

# Interaction q-values are corrected across all three panels at once, exactly as
# the forest plot does it. The text says none passes both q and the sensitivity
# check, and quotes the smallest q (Akkermansia, RS vs Inulin), which fails it.
ic = "DietRS:ConditionIR"
for suffix, why, fmt in (("", "interaction", "q = {:.1e}"), ("_minpres12", "interaction, stricter filter", "({:.1e} under the stricter")):
    tables = load_tables(TAB, f"ANCOMBC2_factorial_Greengenes2{suffix}")
    qj = joint_q(tables)
    ss = lambda panel, t: tables[PANELS[panel][0]].set_index("taxon").loc[t, f"passed_ss_{PANELS[panel][1]}"] in (True, "TRUE", "True")
    hits = [(panel, t, q) for panel, d in qj.items() for t, q in d.items() if q < 0.05 and ss(panel, t)]
    if hits:
        misses.append(f"{why}: text says no interaction passes, tables give {hits}")
    panel, t, q = min(((panel, t, q) for panel, d in qj.items() for t, q in d.items()), key=lambda x: x[2])
    if panel != "RS vs Inulin" or "Akkermansia" not in t or ss(panel, t):
        misses.append(f"{why}: expected Akkermansia (RS vs Inulin, failing the sensitivity check) as the smallest q, got {panel} {t}")
    else:
        expect(fmt.format(q), f"{why} smallest joint q")
ref = load_tables(TAB, "ANCOMBC2_factorial_Greengenes2")["_refInulin"]
akk = ref[ref.taxon.str.contains("Akkermansia")]
expect(f"+{akk[f'lfc_{ic}'].iloc[0] / LN2:.2f} log2 for RS", "interaction estimate")
nominal = sum(int((d[c] < 0.05).sum()) for d, cols in ((main, ["p_DietInulin:ConditionIR", "p_DietRS:ConditionIR"]), (ref, [f"p_{ic}"])) for c in cols)
expect(f"there are {nominal} of them among {3 * len(main)} estimates", "nominal interaction count")
expect(f"{nominal} of {3 * len(main)} estimates nominal where about 5 are expected", "five-ways nominal count")

# ---- untestable-taxa heatmap caption ----
import biom  # noqa: E402
from plot_untestable_heatmap import collapse  # noqa: E402

un = pd.read_csv(os.path.join(TAB, "ANCOMBC2_factorial_Greengenes2_untestable_absent_from_a_cell.tsv"), sep="\t")
sh = pd.read_csv(os.path.join(TAB, "ANCOMBC2_factorial_Greengenes2_untestable_shares.tsv"), sep="\t")
expect(f"The {len(un)} taxa the model cannot test", "heatmap taxon count")
expect(f"{sh.untestable_pct.mean():.2f}% of an average sample's reads", "heatmap mean share")
gm = sh.groupby("Group").untestable_pct.mean()
expect(f"{gm['RS+IR']:.1f}% in RS+IR samples, {gm['Cellulose']:.1f}% in Cellulose and {gm['Cellulose+IR']:.1f}% in Cellulose+IR", "heatmap group shares")
rest = gm[["RS", "Inulin", "Inulin+IR"]]
expect(f"against {rest.min():.1f}–{rest.max():.1f}% in the RS and inulin groups", "heatmap other-group shares")
if gm[["RS", "Inulin", "Inulin+IR"]].max() >= gm[["RS+IR", "Cellulose", "Cellulose+IR"]].min():
    misses.append("heatmap: the three named groups are not the three highest")
fmeta = pd.read_csv(os.path.join(OUT, "metadata_factorial.tsv"), sep="\t", index_col=0, keep_default_na=False)
fc = biom.load_table(os.path.join(os.path.dirname(OUT), "tables", "study-seqs.biom")).to_dataframe(dense=True)[fmeta.index]
ftax = pd.read_csv(os.path.join(os.path.dirname(OUT), "tables", "exported-taxonomy", "Greengenes2_taxonomy.tsv"),
                   sep="\t", index_col=0)["Taxon"]
rel = collapse(fc, ftax) / fc.sum() * 100
row = lambda k: rel.loc[[t for t in rel.index if t.endswith(k)][0]]
am = row("Amulumruptor caecigallinarius"); g = fmeta.Group
if sh.loc[sh.top_taxon_pct.idxmax(), "top_taxon"] != am.name:
    misses.append("heatmap: Amulumruptor is not the largest single share")
ri = am[g == "RS+IR"]
expect(f"{ri.min():.0f}–{ri.max():.0f}% of the reads in every RS+IR sample, at most {am[g != 'RS+IR'].max():.2f}% anywhere else",
       "Amulumruptor shares")
if (ri == 0).any() or am[g.isin(["Inulin", "Inulin+IR"])].gt(0).any():
    misses.append("heatmap: Amulumruptor presence does not match the caption")
ms = row("Mucispirillum schaedleri")
expect(f"reaches {ms.max():.1f}% in one {g[ms.idxmax()]} sample", "Mucispirillum peak")
if ms[g.isin(["Inulin", "Inulin+IR", "RS", "RS+IR"])].gt(0).any():
    misses.append("heatmap: Mucispirillum is not absent from every inulin and RS sample")

# ---- structural zeros: presence table and finding ----
import check_structural_zeros  # noqa: E402
from scipy.stats import fisher_exact  # noqa: E402

sz_miss, szdf = check_structural_zeros.check(html, text)
misses += sz_miss
sizes = fmeta.Group.value_counts()
floor = {c: fisher_exact([[sizes[a], 0], [0, sizes[b]]])[1] for _q, c, a, b in check_structural_zeros.sz.CONTRASTS}
four = [c for _q, c, a, b in check_structural_zeros.sz.CONTRASTS if 4 in (sizes[a], sizes[b])]
if not all(abs(floor[c] - 0.0286) < 1e-3 for c in four) or \
        not all(abs(v - 0.1) < 1e-9 for c, v in floor.items() if c not in four):
    misses.append(f"presence floors differ from the text: {floor}")
for label, n_abs, n in (("RS+IR vs Cellulose+IR", 31, 36), ("RS+IR vs Inulin+IR", 22, 27), ("RS+IR vs RS", 19, 22)):
    d = szdf[szdf.contrast == label]
    if (int((d.absent_from == "RS+IR").sum()), len(d)) != (n_abs, n):
        misses.append(f"direction {label}: tables give {int((d.absent_from == 'RS+IR').sum())} of {len(d)}")
    expect(f"{n_abs} of the {n} against", f"direction {label}")

fig_miss, n_views = check_structural_zeros.check_figures(html, OUT)
misses += fig_miss
if n_views == 0:
    misses.append("no view draws structural zeros")

# ---- presence interaction (presence_interaction.R) ----
def fmtp(p):
    return f"{p:.2g}" if p < 0.01 else (f"{p:.3f}" if p < 0.1 else f"{p:.2g}")


_pi = pd.read_csv(os.path.join(TAB, "ANCOMBC2_factorial_Greengenes2_presence_interaction.tsv"), sep="\t").set_index("term")
_ia, _pd = _pi.loc["Diet:Condition"], _pi.loc["cell"]
expect(f"(p = {fmtp(_ia.p)}, R² {_ia.effect * 100:.0f}%)", "presence interaction p and R2")
expect(f"PERMDISP is borderline (p = {fmtp(_pd.p)})", "presence PERMDISP p")
_nt = len(pd.read_csv(os.path.join(TAB, "ANCOMBC2_factorial_Greengenes2.tsv"), sep="\t"))
_nu = len(pd.read_csv(os.path.join(TAB, "ANCOMBC2_factorial_Greengenes2_untestable_absent_from_a_cell.tsv"), sep="\t"))
expect(f"The plot covers the {_nt} taxa seen in every group. The {_nu} absent", "forest scope")
expect(f"counted over all {_nt + _nu} taxa", "calibration taxon set")
expect(f"Benjamini-Hochberg across the {_nt + _nu}", "Firth taxon set")

# ---- status chips: q < 0.05 without the sensitivity check, or nothing at q < 0.05 ----
def _chip(label):
    m = re.search(rf"(?:<td>|</span>){re.escape(label)}</td>(?:<td[^>]*>[^<]*</td>){{3}}<td><span class=\"chip [a-z]+\">([^<]+)</span>", _da_html)
    return m[1] if m else None
_ss = lambda t, c: t[f"passed_ss_{c}"].astype(str).str.upper().eq("TRUE")
_load = lambda suf: pd.read_csv(os.path.join(TAB, f"ANCOMBC2_factorial_Greengenes2{suf}.tsv"), sep="\t")
_irr = {}
for suf, label in (("_refInulin", "Inulin+IR vs Inulin"), ("_refRS", "RS+IR vs RS")):
    t = _load(suf); c = "ConditionIR"
    f = t[(t[f"q_{c}"] < 0.05) & ~_ss(t, c)]
    _irr[label] = f
    if _chip(label) != ("q &lt; 0.05, not robust" if len(f) else "below threshold"):
        misses.append(f"{label}: chip {_chip(label)!r} but {len(f)} q < 0.05 taxa fail the sensitivity check")
expect(f"On inulin and on RS, {len(_irr['Inulin+IR vs Inulin'])} taxa each reach q", "irradiation near-miss count")
if len(_irr["Inulin+IR vs Inulin"]) != len(_irr["RS+IR vs RS"]):
    misses.append("irradiation near misses differ between inulin and RS, but the text says 2 each")
for label, key, fmt in (("Inulin+IR vs Inulin", "UBA3282 sp003611805", "(up {:.2f} log2 after irradiation)"),
                        ("Inulin+IR vs Inulin", "goldsteinii", "(down {:.2f})"),
                        ("RS+IR vs RS", "Akkermansia", "(up {:.2f} log2 after irradiation, q = "),
                        ("RS+IR vs RS", "UBA3263 sp001689615", "(down {:.2f}), so none")):
    f = _irr[label]; r = f[f.taxon.str.contains(key)]
    if r.empty:
        misses.append(f"{label}: {key} is not a q < 0.05 taxon failing the sensitivity check"); continue
    expect(fmt.format(abs(r["lfc_ConditionIR"].iloc[0]) / LN2), f"{label} {key} estimate")
_ak = _irr["RS+IR vs RS"]; _ak = _ak[_ak.taxon.str.contains("Akkermansia")]
expect(f"q = {_ak.q_ConditionIR.iloc[0]:.1e}) and", "Akkermansia irradiation q on RS (caption)")
expect(f"up {_ak.lfc_ConditionIR.iloc[0] / LN2:.2f} log2 after irradiation on RS (q = {_ak.q_ConditionIR.iloc[0]:.1e})", "Akkermansia irradiation (finding)")
expect(f"rises {_ak.lfc_ConditionIR.iloc[0] / LN2:.2f} log2 on RS", "Akkermansia irradiation (summary)")
_t0 = load_tables(TAB, "ANCOMBC2_factorial_Greengenes2"); _qj = joint_q(_t0)
for panel, d in _qj.items():
    tb = _t0[PANELS[panel][0]].set_index("taxon"); c = PANELS[panel][1]
    n_fail = sum(q < 0.05 and not _ss(tb.loc[[t]], c).iloc[0] for t, q in d.items())
    rows_i = [m for m in re.finditer(rf"<td>{re.escape(panel)}</td>(?:<td[^>]*>[^<]*</td>){{3}}<td><span class=\"chip [a-z]+\">([^<]+)</span>", _da_html)]
    if len(rows_i) != 1 or rows_i[0][1] != ("q &lt; 0.05, not robust" if n_fail else "below threshold"):
        misses.append(f"interaction row {panel}: chip does not match {n_fail} q < 0.05 taxa failing the check")

_akq = min(d[t] for d in _qj.values() for t in d if "Akkermansia" in t)
expect(f"reaches q = {_akq:.4f} but fails the sensitivity check", "five-ways Akkermansia interaction q")

# ---- chance calibration of presence calls (structural_zero_evidence.py) ----
_cal = pd.read_csv(os.path.join(TAB, "ANCOMBC2_factorial_Greengenes2_structural_zeros_calibration.tsv"), sep="\t").set_index(["kind", "contrast"])
def _fp(p):
    return f"{p:.4f}" if p < 0.001 else f"{p:.3f}" if p < 0.1 else f"{p:.2f}"
_i = html.index("Complete separations / chance")
_cal_rows = {c[0]: c[1:] for c in (
    [re.sub(r"<[^>]+>", "", x).strip() for x in re.findall(r"<td[^>]*>(.*?)</td>", r)]
    for r in re.findall(r"<tr[^>]*>(.*?)</tr>", html[_i:html.index("</tbody>", _i)], re.S))}
for (_k, _c), _r in _cal.iterrows():
    if _k != "complete separations":
        continue
    _pa = _cal.loc[("all presence calls", _c)]
    key = "All contrasts" if _c == "all contrasts" else _c
    want = [f"{int(_pa.observed)} / {_pa.expected:.1f}",
            f"{int(_r.observed)} / {_r.expected:.3g}" if _c == "all contrasts" else f"{int(_r.observed)} / {_r.expected:.2g}", _fp(_r.p_perm)]
    if _cal_rows.get(key) != want:
        misses.append(f"calibration row {key}: report {_cal_rows.get(key)}, tables {want}")
_all = _cal.loc[("complete separations", "all contrasts")]
expect(f"Across all contrasts {int(_all.observed)} calls are complete separations, where shuffled labels give {_all.expected:.1f} on average (p = {_fp(_all.p_perm)})", "calibration all")
expect(f"{int(_all.observed)} against {_all.expected:.1f}, p = {_fp(_all.p_perm)}", "summary calibration")
_cs = szdf[szdf.n_present == szdf.n_group]
expect(f"they involve {_cs.taxon.nunique()} taxa", "complete-separation taxa")
for c, other, n_abs in (("RS+IR vs Cellulose+IR", "RS+IR", 18), ("Inulin+IR vs Cellulose+IR", "Inulin+IR", 11)):
    r = _cal.loc[("complete separations", c)]
    if int((_cs[_cs.contrast == c].absent_from == other).sum()) != n_abs:
        misses.append(f"{c}: {n_abs} complete separations missing from {other} in the text, tables differ")
    expect(f"{int(r.observed)} separate {c.split(' vs ')[0]} from Cellulose+IR against {r.expected:.2g} by chance" if c.startswith("RS")
           else f"{int(r.observed)} separate Inulin+IR from Cellulose+IR against {r.expected:.2g} (p = {_fp(r.p_perm)})", f"calibration {c}")
_r = _cal.loc[("complete separations", "RS+IR vs Cellulose+IR")]
expect(f"18 of the {int(_r.observed)} taxa that separate RS+IR from Cellulose+IR completely, against {_r.expected:.2f} expected", "key finding RS calibration")
_r = _cal.loc[("complete separations", "RS+IR vs RS")]
expect(f"{int(_r.observed)} taxa separate RS+IR from unirradiated RS completely against {_r.expected:.1f} by chance (p = {_fp(_r.p_perm)}), "
       f"{int((_cs[_cs.contrast == 'RS+IR vs RS'].absent_from == 'RS+IR').sum())} of them missing after irradiation", "synthesis RS irradiation calibration")
_r = _cal.loc[("complete separations", "Inulin vs Cellulose")]; _r2 = _cal.loc[("complete separations", "RS vs Cellulose")]
expect(f"inulin and cellulose separate on {int(_r.observed)} taxa against {_r.expected:.2g} (p = {_fp(_r.p_perm)})", "calibration Inulin vs Cellulose")
expect(f"{int(_r.observed)} and {int(_r2.observed)} taxa separate completely where chance gives {_r.expected:.2g} and {_r2.expected:.2g} (p = {_fp(_r.p_perm)} and {_fp(_r2.p_perm)})", "synthesis diet calibration")
_cell = lambda k: [int(v) for v in untestable.set_index("taxon").loc[[t for t in untestable.taxon if t.endswith(k)][0]].drop("n_samples")]
for k, want in (("Mucispirillum schaedleri", [3, 0, 0, 4, 0, 0]), ("Bifidobacterium globosum", [0, 3, 2, 0, 4, 0]),
                ("g__Duncaniella;s__", [3, 3, 0, 4, 4, 0]), ("Lactobacillus johnsonii", [0, 3, 0, 4, 1, 0])):
    if _cell(k) != want:
        misses.append(f"{k}: presence per group {_cell(k)}, the text describes {want}")

# ---- Firth per-taxon presence interaction ----
_fr = pd.read_csv(os.path.join(TAB, "ANCOMBC2_factorial_Greengenes2_presence_firth.tsv"), sep="\t")
if not _fr.taxon.iloc[0].endswith("Lactobacillus johnsonii") or (_fr.p_interaction < 0.05).sum() != 1:
    misses.append("Firth: the one raw p < 0.05 taxon is no longer L. johnsonii")
expect(f"absent from every RS sample (p = {_fr.p_interaction.iloc[0]:.4f})", "Firth L. johnsonii p")
expect(f"does not survive correction (q = {_fr.q_interaction.min():.2f})", "Firth smallest q")
expect(f"(p = {_fr.p_interaction.iloc[0]:.4f}, q = {_fr.q_interaction.min():.2f})".replace("0.0098", "0.0098"), "five-ways Firth")

# ---- DA part order ----
_order = [html.find(h) for h in ("How taxa are tested", "Structural zeros: taxa present", "ANCOM-BC2: abundance in the 35 taxa",
                                  "Does irradiation act differently depending on the diet?", ">Both together<", "<h2>Summary</h2>")]
if -1 in _order or _order != sorted(_order):
    misses.append(f"DA parts missing or out of order: {_order}")
_bc = pd.read_csv(os.path.join(TAB, "ANCOMBC2_factorial_Greengenes2_presence_by_cell.tsv"), sep="\t", keep_default_na=False)
_bc = {(r[0], r[1]): r.taxa_present for r in _bc.itertuples(index=False)}
for (_a, _b), _c in [((('Cellulose', 'Control'), ('Cellulose', 'IR')), None), ((('Inulin', 'Control'), ('Inulin', 'IR')), None), ((('RS', 'Control'), ('RS', 'IR')), None)]:
    expect(f"{_bc[_a]:.1f} to {_bc[_b]:.1f}" if _c is None else f"{_bc[_a]:.1f}" , f"presence count {_a}->{_b}")

if misses:
    print("\n".join(misses))
    sys.exit(1)
print(f"report checks pass: {n_fig} figures, all quoted numbers found")
