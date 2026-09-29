#!/usr/bin/env python3
"""Check comparison-Treatment-x-Radiation-factorial/report.html against its tables.

Covers the differential-abundance section and the summary lines that quote it:
every figure the page references exists and is non-empty, the DA table's twelve
rows match da_rows_factorial.py and interaction_rows_factorial.py cell by cell,
and each number in the prose and captions is recomputed from its table. The
diversity sections were audited separately (ANALYSIS_LOG, 2026-09-22) and are
not rechecked here. Exits non-zero listing every miss.

Ported from tisha251224/additional_comparison_scripts/check_report.py.

Usage, inside the qiime2 container:
    python3 check_report_factorial.py
"""
import os
import re
import sys

import numpy as np
import pandas as pd

import interaction_rows_factorial
from da_rows_factorial import markers, res_path, rows
from plot_interaction_forest import coefficient_frame, joint_q, load_panels

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(THIS_DIR), "comparison-Treatment-x-Radiation-factorial")
TAB = os.path.join(OUT, "tables")
STEM = "ANCOMBC2_factorial_Greengenes2_Day17_minpres6"

html = open(os.path.join(OUT, "report.html")).read()
text = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", html))
misses = []


def expect(value, why):
    if value not in text:
        misses.append(f"{why}: '{value}' not in report")


def cells(row_html):
    return [re.sub(r"<[^>]+>", "", c).strip() for c in re.findall(r"<td[^>]*>(.*?)</td>", row_html)]


def fmt_row(n, up, down, q):
    u = f"+{up:.2f}" if up is not None else "—"
    d = f"−{down:.2f}" if down is not None else "—"
    return (str(n), "—" if n == 0 else f"{u} / {d}", "—" if q is None else f"{q:.1e}")


# ---- figures ----
for src in set(re.findall(r'src="([^"]+)"', html)):
    path = os.path.join(OUT, src)
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        misses.append(f"figure missing or empty: {src}")
n_fig = len(set(re.findall(r'src="(figures/[^"]+\.png)"', html)))
expect(f"{n_fig} figures, PNG and SVG", "figure count")

# ---- tested set ----
main = pd.read_csv(os.path.join(TAB, f"{STEM}.tsv"), sep="\t")
strict = pd.read_csv(os.path.join(TAB, "ANCOMBC2_factorial_Greengenes2_Day17_minpres12.tsv"), sep="\t")
named = main.taxon.map(lambda t: bool(re.match(r"s__.", t.split(";")[-1].strip())))
expect(f"{len(main)} of them are seen in at least one mouse of every one of the six groups", "tested count")
expect(f"{named.sum()} carry a species call, and the other {(~named).sum()}", "species split")
expect(f"and here there are {len(main)}", "taxa vs the 50 recommendation")
expect(f"least 12 of the 24 mice ({len(strict)} tested)", "stricter filter taxa")
_n_un = len(pd.read_csv(os.path.join(TAB, f"{STEM}_untestable_absent_from_a_cell.tsv"), sep="\t"))
expect(f"splits its power across the {len(main)} taxa seen in every group, and the presence tests cover the {_n_un} it cannot use", "summary taxa count")

# ---- untestable taxa and heatmap ----
un = pd.read_csv(os.path.join(TAB, f"{STEM}_untestable_absent_from_a_cell.tsv"), sep="\t")
sh = pd.read_csv(os.path.join(TAB, f"{STEM}_untestable_shares.tsv"), sep="\t")
absent = un[[c for c in un.columns if "(n=" in c]].eq(0)
n_ifa = int(absent["IFA (n=4)"].sum())
expect(f"{len(un)} taxa are seen in at least 6 mice but are missing entirely", "untestable count")
expect(f"The other {len(un)} are missing from at least one whole group", "structural zero count")
expect(f"{n_ifa} of the {len(un)} are absent from every unirradiated IFA mouse", "absent from IFA")
expect(f"The {len(un)} taxa that could not be tested", "heatmap taxon count")
expect(f"{int((absent.sum(axis=1) == 1).sum())} are missing from one group only", "one-group count")
gm = sh.groupby("Group").untestable_pct.mean()
expect(f"{sh.untestable_pct.mean():.2f}% of an average mouse's reads, but {gm['IFA']:.2f}% in unirradiated IFA mice, "
       f"where {n_ifa} of the {len(un)} are absent", "heatmap shares")
if gm.idxmin() != "IFA":
    misses.append("heatmap: IFA is not the group with the lowest untestable share")
tt = lambda k: un.set_index("taxon").loc[[t for t in un.taxon if t.endswith(k)][0]]
top = sh.loc[sh.top_taxon_pct.idxmax()]
ac = tt("acidifaciens")
if "acidifaciens" not in top.top_taxon or top.Group != "Control":
    misses.append(f"heatmap: largest single share is {top.top_taxon} in {top.Group}")
expect(f"acidifaciens reaches {top.top_taxon_pct:.1f}% of a Control mouse's reads and is absent from "
       f"{['no', 'one', 'two', 'three', 'four', 'five', 'six'][int(absent.loc[un.taxon.str.endswith('acidifaciens')].sum(axis=1).iloc[0])]} of the six groups",
       "B. acidifaciens")
al_row = un.taxon.str.endswith("Alloprevotella sp002933955")
al_abs = [c.split(" (n=")[0] for c in absent.columns if absent.loc[al_row, c].iloc[0]]
al_max = sh.top_taxon_pct[sh.top_taxon == un.taxon[al_row].iloc[0]].max()
if al_abs != ["Control", "IFA"]:
    misses.append(f"heatmap: Alloprevotella is absent from {al_abs}")
expect(f"sp002933955 reaches {al_max:.1f}% and is absent from both Control and IFA", "Alloprevotella")
ca = tt("caccae")
if ca["Control (n=4)"] != 0 or ca["n_samples"] != 20:
    misses.append("heatmap: B. caccae is not in exactly the 20 non-Control mice")

# ---- DA table ----
table = {}
question = None
i = html.index("Prevalence ≥ 12")
body = html[html.index("<tbody>", i):html.index("</tbody>", i)]
for r in re.findall(r"<tr[^>]*>(.*?)</tr>", body, re.S):
    c = cells(r)
    if len(c) == 7:
        question, c = c[0], c[1:]
    table[(question, c[0])] = c
for r in rows():
    got = table.get((r["question"], r["label"]))
    want = fmt_row(r["n"], r["max_up"], r["max_down"], r["min_q"]) + (str(r["n_strict"]),)
    if got is None:
        misses.append(f"DA row {r['label']}: not found")
    elif tuple(got[1:5]) != want:
        misses.append(f"DA row {r['label']}: report {tuple(got[1:5])}, tables {want}")
strict_rows = {r["label"]: r for r in interaction_rows_factorial.rows("_Day17_minpres12")}
for r in interaction_rows_factorial.rows():
    got = table.get(("Does radiation's effect differ between treatments?", r["label"]))
    want = fmt_row(r["n"], r["max_up"], r["max_down"], r["min_q"]) + (str(strict_rows[r["label"]]["n"]),)
    if got is None or tuple(got[1:5]) != want:
        misses.append(f"interaction row {r['label']}: report {got and tuple(got[1:5])}, tables {want}")

# ---- prose on the DA rows ----
R = {r["label"]: r for r in rows()}
ifa, ino, ifa_ino = R["IFA vs Control"], R["Inosine vs Control"], R["IFA vs Inosine"]
expect(f"−{ino['max_down']:.2f} for Inosine vs Control means a taxon 2 {ino['max_down']:.2f} , "
       f"about {2 ** ino['max_down']:.0f}-fold", "fold-change example")
if ifa["n"] != 0 or ifa_ino["n"] != 0 or ifa_ino["n_strict"] != 0:
    misses.append("IFA finding says IFA vs Control and IFA vs Inosine have no hit at the default filter "
                  "and IFA vs Inosine none under the stricter one")
expect(f"Under the stricter prevalence filter {ifa['n_strict']} pass, all in", "IFA stricter-filter count")
strict_up = [m[0] for m in markers(res_path("_minpres12", "_pairwise_TreatmentIFA"))]
if not all("f__Muribaculaceae" in x for x in strict_up) or \
        any(m[1] != "IFA" for m in markers(res_path("_minpres12", "_pairwise_TreatmentIFA"))):
    misses.append("not every stricter-filter IFA vs Control hit is a Muribaculaceae taxon higher under IFA")
ir = R["IFA+IR vs IR"]
if ir["n"] != 1 or ir["down_taxa"] != ["s__Lawsonibacter_sp000492175"]:
    misses.append(f"IFA+IR vs IR: hits are {ir['up_taxa'] + ir['down_taxa']}, prose says Lawsonibacter alone")
expect(f"(down {ir['max_down']:.2f}, q = {ir['min_q']:.1e})", "IFA+IR vs IR estimate and q")
expect(f"Inosine differs from Control in {ino['n']} taxon", "inosine count")
expect(f"({ino['n_strict']} under the stricter filter)", "inosine stricter count")
if ino["down_taxa"] != ["s__Alistipes_A_871400_sp"]:
    misses.append(f"Inosine vs Control: hit is {ino['down_taxa']}, prose says Alistipes_A")
expect(f"(down {ino['max_down']:.2f}, q = {ino['min_q']:.1e})", "inosine estimate and q")
# Counts of taxa with q < 0.05 that fail the sensitivity check, as the captions give them.
FIT = os.path.join(TAB, STEM)
for suf, coef, label in (("_pairwise_pairs", "TreatmentIFA", "IFA vs Control"),
                         ("_pairwise_pairs", "TreatmentInosine", "Inosine vs Control"),
                         ("_pairwise_pairs", "TreatmentInosine_TreatmentIFA", "IFA vs Inosine"),
                         ("_radrefIR_pairwise_pairs", "TreatmentInosine", "Inosine+IR vs IR"),
                         ("_radrefIR_pairwise_pairs", "TreatmentInosine_TreatmentIFA", "IFA+IR vs Inosine+IR"),
                         ("", "RadiationIR", "within Control"), ("_refIFA", "RadiationIR", "within IFA")):
    t = pd.read_csv(f"{FIT}{suf}.tsv", sep="\t")
    n_fail = int(((t[f"q_{coef}"] < 0.05) & ~t[f"passed_ss_{coef}"].astype(str).str.upper().eq("TRUE")).sum())
    expect(f"{n_fail} reach q &lt; 0.05 but fail the sensitivity check", f"{label}: failing count")

# ---- interaction caption ----
tables, panels = load_panels(TAB, STEM)
qj = joint_q(tables, panels)
irows = {r["label"]: r for r in interaction_rows_factorial.rows()}
if any(r["n"] for r in irows.values()):
    misses.append("interaction: text says no coefficient passes; the tables give hits")
nom = [irows[k]["nominal"] for k in ("IFA vs Control", "Inosine vs Control", "IFA vs Inosine")]
words = ["no", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten", "eleven", "twelve"]
expect(f"With {len(main)} taxa per panel, about {words[round(0.05 * len(main))]} such hits per panel", "chance count")
expect(f"The IFA panel has {words[nom[0]]}", "IFA panel nominal")
expect(f"the Inosine panel has {words[nom[1]]}", "Inosine panel nominal")
expect(f"the IFA-vs-Inosine panel {words[nom[2]]}", "third panel nominal")
expect(f"{sum(nom)} of {3 * len(main)} coefficients nominal where about {round(0.05 * 3 * len(main))} are expected", "summary nominal count")
failing = []
for panel, (fit, coef, sign) in panels.items():
    df = tables[fit].set_index("taxon")
    failing += [(panel, t, q) for t, q in qj[panel].items()
                if q < 0.05 and str(df.loc[t, f"passed_ss_{coef}"]).upper() != "TRUE"]
expect(f"{words[len(failing)].capitalize()} reach q &lt; 0.05 but fail the check", "failing interaction count")
for panel, t, q in failing:
    expect(f"q = {q:.1e}", f"interaction {panel} {t.split(';')[-1]} q")
g = main.taxon.map(lambda t: [p.strip() for p in t.split(";") if p.strip() and not p.strip().endswith("__")])
expect(f"limited to the {int(g.map(lambda ps: bool(ps) and not ps[-1].startswith('s__')).sum())} taxa named", "genus forest count")

# ---- structural zeros: presence table and finding ----
import check_structural_zeros  # noqa: E402
from scipy.stats import fisher_exact  # noqa: E402

check_structural_zeros.sz.configure("factorial")
sz_miss, szdf = check_structural_zeros.check(html, text)
misses += sz_miss
expect(f"gives p = {fisher_exact([[4, 0], [0, 4]])[1]:.3f} at best", "design floor")
for label, grp, n_abs, n in (("IFA vs Control", "IFA", 17, 18), ("IFA vs Inosine", "IFA", 20, 23), ("IFA", "IFA", 17, 20)):
    d = szdf[szdf.contrast == label]
    if (int((d.absent_from == grp).sum()), len(d)) != (n_abs, n):
        misses.append(f"direction {label}: tables give {int((d.absent_from == grp).sum())} of {len(d)}")
    expect(f"{n_abs} of the {n}", f"direction {label}")

fig_miss, n_views = check_structural_zeros.check_figures(html, OUT)
misses += fig_miss
if n_views == 0:
    misses.append("no view draws structural zeros")

# ---- presence interaction (presence_interaction.R) ----
def fmtp(p):
    return f"{p:.2g}" if p < 0.01 else (f"{p:.3f}" if p < 0.1 else f"{p:.2g}")


_pi = pd.read_csv(os.path.join(TAB, "ANCOMBC2_factorial_Greengenes2_Day17_minpres6_presence_interaction.tsv"), sep="\t").set_index("term")
_ia, _pd = _pi.loc["Treatment:Radiation"], _pi.loc["cell"]
expect(f"(p = {fmtp(_ia.p)}, R² {_ia.effect * 100:.0f}%)", "presence interaction p and R2")
expect(f"PERMDISP p = {fmtp(_pd.p)}", "presence PERMDISP p")
_nt = len(pd.read_csv(os.path.join(TAB, "ANCOMBC2_factorial_Greengenes2_Day17_minpres6.tsv"), sep="\t"))
_nu = len(pd.read_csv(os.path.join(TAB, "ANCOMBC2_factorial_Greengenes2_Day17_minpres6_untestable_absent_from_a_cell.tsv"), sep="\t"))
expect(f"The plot covers the {_nt} taxa seen in every group. The {_nu} absent", "forest scope")
expect(f"{_nt + _nu} taxa are seen in at least 6 of the 24 mice", "presence taxon set")
_bc = pd.read_csv(os.path.join(TAB, "ANCOMBC2_factorial_Greengenes2_Day17_minpres6_presence_by_cell.tsv"), sep="\t", keep_default_na=False)
_bc = {(r[0], r[1]): r.taxa_present for r in _bc.itertuples(index=False)}
for (_a, _b), _c in [((('IFA', 'None'), ('IFA', 'None')), 1), ((('Inosine', 'None'), ('Inosine', 'None')), 1), ((('Control', 'None'), ('Control', 'None')), 1), ((('IFA', 'IR'), ('IFA', 'IR')), 1), ((('Control', 'IR'), ('Control', 'IR')), 1), ((('Inosine', 'IR'), ('Inosine', 'IR')), 1)]:
    expect(f"{_bc[_a]:.1f} to {_bc[_b]:.1f}" if _c is None else f"{_bc[_a]:.1f}" , f"presence count {_a}->{_b}")

# ---- structural-zero evidence (structural_zero_evidence.py) ----
_S = os.path.join(TAB, "ANCOMBC2_factorial_Greengenes2_Day17_minpres6")
_cal = pd.read_csv(_S + "_structural_zeros_calibration.tsv", sep="\t").set_index("contrast")
_dep = pd.read_csv(_S + "_structural_zeros_depth.tsv", sep="\t")
_t = _cal.loc["all contrasts"]
expect(f"{_t.observed} calls in {_dep.taxon.nunique()} taxa show complete separation, where shuffled labels give "
       f"{_t.expected:.1f} on average (permutation p = {_t.p_perm:.4f})", "SZ calibration")
if (_dep.q_depth >= 0.05).any():
    misses.append("text says no complete-separation absence is explained by depth")
expect(f"{_dep.taxon.nunique()} separate groups completely ({_t.observed} calls against {_t.expected:.1f} expected", "SZ summary row")
_rad = _dep.contrast.isin(["Control", "IFA", "Inosine"]).sum()
if (round(3 * _rad / len(_dep)), round(3 * (len(_dep) - _rad) / len(_dep))) != (1, 2):
    misses.append(f"text says a third of the calls are irradiation within a treatment; tables give {_rad} of {len(_dep)}")
if _dep.contrast.value_counts()[["Control", "IFA", "Inosine"]].idxmax() != "Control":
    misses.append("text says irradiation separates the most in Control mice")
_un = pd.read_csv(_S + "_untestable_absent_from_a_cell.tsv", sep="\t").set_index("taxon")
_cc = [c for c in _un.columns if "(n=" in c]
_rows = re.findall(r"<tr><td>(.*?)</td><td>([^<]*)</td><td>([^<]*)</td><td class=\"n\">([0-9.]+)%</td></tr>",
                   html[html.index("In every mouse of"):])
if len(_rows) != _dep.taxon.nunique():
    misses.append(f"SZ table has {len(_rows)} rows, tables give {_dep.taxon.nunique()} taxa")
_want = []
for t in _dep.taxon.unique():
    u = _un.loc[t, _cc]
    _want.append((", ".join(c.split(" (")[0] for c in _cc if u[c] == 4), ", ".join(c.split(" (")[0] for c in _cc if u[c] == 0),
                  f"{_dep.mean_share_pct[_dep.taxon == t].max():.2f}"))
if sorted(_want) != sorted((r[1], r[2], r[3]) for r in _rows):
    misses.append("SZ table rows (groups and shares) do not match the tables")
for _k in ("Alloprevotella sp002933955", "Desulfovibrio sp011039135"):
    _u = _un.loc[[t for t in _un.index if t.endswith(_k)][0]]
    if not (_u["IR (n=4)"] == 4 and _u["Control (n=4)"] == 0):
        misses.append(f"{_k}: not in every IR mouse and no Control mouse")
_ac = _un.loc[[t for t in _un.index if t.endswith("acidifaciens")][0], _cc]
if [c.split(" (")[0] for c in _cc if _ac[c] == 0] != ["IFA", "Inosine", "IR", "Inosine+IR"] or _ac["Control (n=4)"] != 4:
    misses.append("B. acidifaciens is not in every Control mouse and absent from IR, IFA, Inosine and Inosine+IR")

# ---- restructured DA section: Firth per taxon, and the synthesis claims ----
_sz = pd.read_csv(os.path.join(TAB, f"{STEM}_structural_zeros.tsv"), sep="\t")
expect(f"all {len(_sz)} such tests", "SZ test count")
_fr = pd.read_csv(os.path.join(TAB, f"{STEM}_presence_firth.tsv"), sep="\t")
expect(f"{int((_fr.p_interaction < 0.05).sum())} of the {len(_fr)} have raw p &lt; 0.05 and one passes correction", "Firth counts")
if int((_fr.q_interaction < 0.05).sum()) != 1 or not _fr.taxon.iloc[0].endswith("g__Duncaniella;s__"):
    misses.append("Firth: text says exactly one taxon passes, the genus-level Duncaniella")
expect(f"Duncaniella (q = {_fr.q_interaction.iloc[0]:.2g})", "Firth Duncaniella q")
expect(f"(q = {_fr.q_interaction.iloc[1]:.2f}), is in every", "Firth second taxon q")
if "acidifaciens" not in _fr.taxon.iloc[1]:
    misses.append("Firth: second-strongest taxon is not B. acidifaciens")
_c = lambda k: _un.loc[[t for t in _un.index if t.endswith(k)][0], _cc].rename(lambda c: c.split(" (")[0])
_dn = _c("g__Duncaniella;s__")
if tuple(_dn[["Control", "IR", "IFA", "IFA+IR", "Inosine", "Inosine+IR"]]) != (0, 4, 0, 4, 4, 0):
    misses.append(f"Duncaniella presence is not 0->4 in Control and IFA, 4->0 in Inosine: {dict(_dn)}")
_ba = _c("acidifaciens")
if tuple(_ba[["Control", "IR", "IFA", "IFA+IR", "Inosine", "Inosine+IR"]]) != (4, 0, 0, 3, 0, 0):
    misses.append(f"B. acidifaciens presence differs from the text: {dict(_ba)}")
_sh = lambda t: t.split(";")[-1] if not t.endswith("__") else t.split(";")[-2]
_ino = {(_sh(r.taxon), r.present_in.replace("Inosine", "X")) for r in _dep[_dep.contrast == "Inosine vs Control"].itertuples()}
_irc = {(_sh(r.taxon), r.present_in.replace("IR", "X")) for r in _dep[_dep.contrast == "Control"].itertuples()}
if len(_ino) != 4 or not _ino <= _irc:
    misses.append("text says the four Inosine-vs-Control complete separations are all IR-vs-Control ones, same direction")
_ifa_c = _sz[_sz.contrast == "IFA vs Control"]
expect(f"{int((_ifa_c.absent_from == 'IFA').sum())} of the {len(_ifa_c)} presence calls against Control are taxa missing from IFA", "IFA presence calls")
expect(f"IFA mice carry {_bc[('IFA', 'None')]:.1f} of the 79 taxa against {_bc[('Control', 'None')]:.1f} in Control", "IFA presence count")
for _h in ("How taxa are tested", "Structural zeros: taxa present in some groups", "ANCOM-BC2: abundance in the 46 taxa", "Both together"):
    if _h not in html:
        misses.append(f"DA part heading missing: {_h}")
_order = [html.find(_h) for _h in ("How taxa are tested", "Structural zeros: taxa present", "ANCOM-BC2: abundance in", "Both together", "<h2>Summary</h2>")]
if _order != sorted(_order):
    misses.append("DA parts are out of order")

if misses:
    print("\n".join(misses))
    sys.exit(1)
print(f"report checks pass: {n_fig} figures, all quoted numbers found")
