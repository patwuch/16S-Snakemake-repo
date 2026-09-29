#!/usr/bin/env python3
"""Check comparison-Diet-x-Radiation-factorial/report.html against its tables.

Every figure the page references must exist and be non-empty, and each number
below is recomputed from its table and must appear in the report as formatted,
in its own table row where it has one. Exits non-zero listing every miss.

Ported from tisha251224/additional_comparison_scripts/check_report.py.

Usage, inside the qiime2 container:
    python3 check_report.py
"""
import math
import os
import re
import sys

import biom
import numpy as np
from scipy.stats import fisher_exact
import pandas as pd

from da_rows import SWEEP, rows
from plot_interaction_forest import PANELS, joint_q, load_tables
from plot_untestable_heatmap import collapse

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(THIS_DIR), "comparison-Diet-x-Radiation-factorial")
TAB = os.path.join(OUT, "tables")
LN2 = np.log(2)
PRE = "ANCOMBC2_factorial_Greengenes2"

html = open(os.path.join(OUT, "report.html")).read()
text = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", html))
misses = []


def expect(value, why):
    if value not in text:
        misses.append(f"{why}: '{value}' not in report")


def sig2(v):
    """Two significant figures, trailing zero kept: 0.0302 -> '0.030'."""
    return f"{v:.{max(0, 1 - math.floor(math.log10(abs(v))))}f}"


def cells(row_html):
    return [re.sub(r"<[^>]+>", "", c).strip() for c in re.findall(r"<td[^>]*>(.*?)</td>", row_html)]


def table_rows(after):
    """The <tr> rows of the first table following the heading text `after`."""
    i = html.index(after)
    body = html[html.index("<tbody>", i):html.index("</tbody>", i)]
    return [cells(r) for r in re.findall(r"<tr[^>]*>(.*?)</tr>", body, re.S)]


# ---- figures ----
for src in set(re.findall(r'src="([^"]+)"', html)):
    path = os.path.join(OUT, src)
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        misses.append(f"figure missing or empty: {src}")
n_fig = len(set(re.findall(r'src="(figures/[^"]+\.png)"', html)))
expect(f"{n_fig} figures, PNG and SVG", "figure count")

# ---- design ----
meta = pd.read_csv(os.path.join(OUT, "metadata_factorial.tsv"), sep="\t")
expect(f"Samples {len(meta)}", "sample count")
per_cell = meta.groupby(["Diet", "Radiation"]).size()
if per_cell.nunique() != 1:
    misses.append(f"design is not balanced: {per_cell.to_dict()}")
expect(f"n = {per_cell.iloc[0]} mice", "per-cell n")

# ---- diversity ----
div = pd.read_csv(os.path.join(TAB, "diversity_factorial.tsv"), sep="\t")
alpha = div[div.analysis == "alpha"].set_index(["metric", "term"])
METRIC = {"Faith PD": "faith_pd", "Chao1": "chao1", "Shannon": "shannon"}
TERM = {"Diet × Radiation": "Diet:Radiation", "Diet": "Diet"}
for r in table_rows("<h2>Alpha diversity</h2>"):
    a = alpha.loc[(METRIC[r[0]], TERM[r[1]])]
    if (r[3], r[4]) != (sig2(a.p_raw), sig2(a.p_BH)):
        misses.append(f"alpha {r[0]} {r[1]}: report ({r[3]}, {r[4]}), table ({sig2(a.p_raw)}, {sig2(a.p_BH)})")
expect(f"smallest q is {sig2(alpha.p_BH.min())}", "alpha smallest q")
# The table must list the lowest q-values, as the text says.
shown = sorted(float(r[4]) for r in table_rows("<h2>Alpha diversity</h2>"))
if shown != sorted(float(sig2(q)) for q in alpha.p_BH.nsmallest(len(shown))):
    misses.append("alpha table does not list the lowest q-values")
sh = alpha.groupby(level="metric").resid_normal_p.first()
expect(f"Shapiro p = {sh['simpson']:.1g} and {sh['evenness']:.2g}", "Shapiro p")

# Group means quoted in the direction column.
av = pd.read_csv(os.path.join(TAB, "diversity_alpha_values.tsv"), sep="\t", keep_default_na=False, na_values=[""])
cm = av.groupby(["metric", "Diet", "Radiation"]).value.mean()
fmt = {"faith_pd": lambda v: f"{v:,.0f}", "chao1": lambda v: f"{v:.0f}", "shannon": lambda v: f"{v:.2f}"}
for m, f in fmt.items():
    for d in ("Control", "SCFAs"):
        expect(f"{f(cm[m, d, 'None'])} → {f(cm[m, d, 'IR'])}", f"alpha {m} {d} IR shift")
for m in ("faith_pd", "chao1"):
    dm = cm[m].groupby(level="Diet").mean()
    expect(f"({fmt[m](dm['Control'])} → {fmt[m](dm['SCFAs'])})", f"alpha {m} diet means")

BETA = {"Jaccard": "jaccard", "Bray-Curtis": "bray_curtis",
        "Unweighted UniFrac · Greengenes2": "unweighted_unifrac",
        "Weighted UniFrac · Greengenes2": "weighted_unifrac"}
perm = div[div.analysis == "beta_permanova"].set_index(["term", "metric"])
term = None
for r in table_rows("PERMANOVA: do the groups differ in location?"):
    if len(r) == 6:
        term, r = r[0], r[1:]
    t = {"Interaction": "Diet:Radiation"}.get(term, term)
    p = perm.loc[(t, BETA[r[0]])]
    want = (sig2(p.p_raw), sig2(p.p_BH), f"{p.effect * 100:.0f}%")
    if (r[1], r[2], r[4]) != want:
        misses.append(f"PERMANOVA {term} {r[0]}: report {(r[1], r[2], r[4])}, table {want}")
    passed = "not" not in r[3]
    if passed != (p.p_BH < 0.05):
        misses.append(f"PERMANOVA {term} {r[0]}: result chip disagrees with q")
eff = lambda t: perm.xs(t, level="term").effect * 100
for t, why in (("Diet", "Diet"), ("Diet:Radiation", "interaction"), ("Residual", "residual")):
    expect(f"{eff(t).min():.0f}–{eff(t).max():.0f}%", f"{why} R2 range")
expect(f"irradiation {eff('Radiation').min():.0f}%", "radiation R2")

DISP = {"Cell (all four groups)": ("beta_dispersion", "Cell"), "Diet": ("beta_dispersion_factorial", "Diet"),
        "Radiation": ("beta_dispersion_factorial", "Radiation"),
        "Interaction": ("beta_dispersion_factorial", "Diet:Radiation")}
disp = div.set_index(["analysis", "term", "metric"])
for r in table_rows("PERMDISP: do the groups differ in spread?"):
    if len(r) == 4:
        term, r = r[0], r[1:]
    q = disp.loc[(*DISP[term], BETA[r[0]])].p_BH
    if r[1] != sig2(q):
        misses.append(f"PERMDISP {term} {r[0]}: report {r[1]}, table {sig2(q)}")
    if q < 0.05:
        misses.append(f"PERMDISP {term} {r[0]}: q < 0.05 but shown as equal spread")

# ---- differential abundance ----
main = pd.read_csv(os.path.join(TAB, f"{PRE}.tsv"), sep="\t")
untestable = pd.read_csv(os.path.join(TAB, f"{PRE}_untestable_absent_from_a_cell.tsv"), sep="\t")
named = main.taxon.map(lambda t: bool(re.match(r"s__.", t.split(";")[-1].strip())))
expect(f"{len(main)} of them are seen in at least one mouse of every one of the four groups", "tested count")
expect(f"{named.sum()} carry a species name, and the other {(~named).sum()}", "species split")
expect(f"The other {len(untestable)} are missing from at least one whole group", "untestable count")
sir = untestable.columns[untestable.columns.str.startswith("SCFAs+IR")][0]
# Untestable-taxa heatmap caption.
cells_un = [c for c in untestable.columns if "(n=" in c]
absent = untestable[cells_un].eq(0)
one = absent.sum(axis=1) == 1
for g in ("SCFAs+IR", "Control", "SCFAs", "IR"):
    col = next(c for c in cells_un if c.startswith(f"{g} (n="))
    k = int((one & absent[col]).sum())
    phrase = f"{k} are missing only from {g}" if g == "SCFAs+IR" else f"{k} only from {g}"
    expect(phrase, f"untestable missing only from {g}")
expect(f"the other {int((absent.sum(axis=1) == 2).sum())} are missing from two groups", "untestable two-group count")
if (absent.sum(axis=1) > 2).any():
    misses.append("an untestable taxon is missing from more than two groups; caption does not cover it")
shares = pd.read_csv(os.path.join(TAB, f"{PRE}_untestable_shares.tsv"), sep="\t")
expect(f"{shares.untestable_pct.mean():.2f}% of an average mouse's reads", "untestable mean share")
rod = untestable.taxon[untestable.taxon.str.contains("rodentium")].iloc[0]
counts = biom.load_table(os.path.join(os.path.dirname(OUT), "tables", "study-seqs.biom")).to_dataframe(dense=True)
tax = pd.read_csv(os.path.join(os.path.dirname(OUT), "tables", "exported-taxonomy", "Greengenes2_taxonomy.tsv"),
                  sep="\t", index_col=0)["Taxon"]
has_rod = collapse(counts, tax).loc[rod] > 0
prefixes = sorted({m.split(".")[0] for m in has_rod.index[has_rod]})
whole = all(has_rod[[m for m in has_rod.index if m.split(".")[0] == p]].all() for p in prefixes)
if prefixes != ["N4", "N5", "N6"] or not whole:
    misses.append(f"rodentium presence is not exactly every N4, N5 and N6 mouse: {sorted(has_rod.index[has_rod])}")
expect(f"up to {shares.top_taxon_pct[shares.top_taxon == rod].max():.1f}% of a mouse's reads, is in "
       f"{int(untestable.set_index('taxon').loc[rod, 'n_samples'])} mice", "rodentium share and prevalence")
fl = shares.loc[shares.top_taxon_pct.idxmax()]
expect(f"reaches {fl.top_taxon_pct:.1f}% in {fl.sampleid}", "flavefaciens peak")
if "flavefaciens" not in fl.top_taxon:
    misses.append(f"largest single untestable share is {fl.top_taxon}, not flavefaciens_B")
for n in SWEEP:
    t = pd.read_csv(os.path.join(TAB, f"{PRE}_minpres{n}.tsv"), sep="\t")
    expect(f"{len(t)}" + (" taxa tested" if n == SWEEP[-1] else ""), f"taxa tested at >= {n}")
expect(", ".join(str(len(pd.read_csv(os.path.join(TAB, f"{PRE}_minpres{n}.tsv"), sep="\t"))) for n in SWEEP[:2]),
       "sweep taxa counts")

# DA table rows keyed by (question, contrast): the unirradiated diet row and the
# interaction row share the label "SCFAs vs Control".
da, question = {}, None
for r in table_rows("Prevalence ≥ 8 / 12 / 16"):
    if len(r) == 7:
        question, r = r[0], r[1:]
    da[(question, r[0])] = r
for r in rows():
    got = da.get((r["question"], r["label"]))
    if got is None:
        misses.append(f"DA row {r['label']}: not found")
        continue
    up = f"+{r['max_up']:.2f}" if r["max_up"] is not None else "—"
    down = f"−{r['max_down']:.2f}" if r["max_down"] is not None else "—"
    ud = "—" if r["n"] == 0 else f"{up} / {down}"
    want = (str(r["n"]), ud, "—" if r["min_q"] is None else f"{r['min_q']:.1e}",
            " / ".join(str(r["sweep"][n]) for n in SWEEP))
    if tuple(got[1:5]) != want:
        misses.append(f"DA row {r['label']}: report {tuple(got[1:5])}, tables {want}")

# The interaction row: one panel, q by one BH pass (joint_q), hit = q < 0.05 and passed_ss.
ic = PANELS["SCFAs vs Control"][1]


def interaction_hits(suffix):
    tables = load_tables(TAB, f"{PRE}{suffix}")
    df = tables[""].set_index("taxon")
    qj = joint_q(tables)["SCFAs vs Control"]
    ok = lambda t: str(df.loc[t, f"passed_ss_{ic}"]).upper() == "TRUE"
    return df, qj, [t for t, q in qj.items() if q < 0.05 and ok(t)]


df, qj, hits = interaction_hits("")
if hits:
    misses.append(f"interaction: text says none passes, tables give {hits}")
row = da[("Does radiation's effect differ between diets?", "SCFAs vs Control")]
want = ("0", "—", "—", " / ".join(str(len(interaction_hits(f'_minpres{n}')[2])) for n in SWEEP))
if tuple(row[1:5]) != want:
    misses.append(f"DA interaction row: report {tuple(row[1:5])}, tables {want}")
dun = [t for t in df.index if t.rstrip("; _s").endswith("g__Duncaniella")][0]
if str(df.loc[dun, f"passed_ss_{ic}"]).upper() == "TRUE" or qj[dun] >= 0.05:
    misses.append("interaction: genus Duncaniella is no longer a q-passing taxon that fails the sensitivity check")
expect(f"+{df.loc[dun, f'lfc_{ic}'] / LN2:.2f} log2 (q = {qj[dun]:.3f})", "Duncaniella interaction estimate and q")
_ix = lambda k: [t for t in df.index if t.endswith(k)][0]
for _k, _txt in (("Acetatifactor sp011959105", "sp011959105 at +"), (";g__Dysosmobacter", "Dysosmobacter at −")):
    _t = _ix(_k)
    if not qj[_t] < 0.05 or str(df.loc[_t, f"passed_ss_{ic}"]).upper() == "TRUE":
        misses.append(f"interaction: {_k} is no longer q < 0.05 and failing the sensitivity check")
    expect(f"{_txt}{abs(df.loc[_t, f'lfc_{ic}']) / LN2:.2f} log2", f"{_k} interaction estimate")
# Same coefficient, up to sign, in the two refits.
for suf in ("_radrefIR", "_refSCFAs"):
    o = pd.read_csv(os.path.join(TAB, f"{PRE}{suf}.tsv"), sep="\t").set_index("taxon")
    c = [k for k in o.columns if k.startswith("lfc_") and ":" in k][0]
    if not np.allclose(o[c].reindex(df.index).abs(), df[f"lfc_{ic}"].abs(), equal_nan=True):
        misses.append(f"interaction differs in {suf} fit")
# Per-row count of taxa with q < 0.05 that fail the sensitivity check, as captions and finding give them.
_fail = {}
for suf, coef, label in (("", "DietSCFAs", "SCFAs vs Control"), ("_radrefIR", "DietSCFAs", "SCFAs+IR vs IR"),
                         ("", "RadiationIR", "IR vs Control"), ("_refSCFAs", "RadiationIR", "SCFAs+IR vs SCFAs")):
    t = pd.read_csv(os.path.join(TAB, f"{PRE}{suf}.tsv"), sep="\t")
    n_fail = int(((t[f"q_{coef}"] < 0.05) & ~t[f"passed_ss_{coef}"].astype(str).str.upper().eq("TRUE")).sum())
    expect(f"{n_fail} taxa reach q &lt; 0.05 but fail the sensitivity check", f"{label}: q-passing taxa failing the sensitivity check (caption)")
    _fail[label] = n_fail
p = pd.to_numeric(df[f"p_{ic}"], errors="coerce")
nominal = int((p < 0.05).sum())
expect(f"Of the {len(df)} taxa, {nominal} have raw p &lt; 0.05", "nominal interaction count")
expect(f"about {round(0.05 * len(df))} would be expected", "expected by chance")
q_ok = df.index[(pd.Series(qj) < 0.05).reindex(df.index).fillna(False)]
_fail_ix = len(q_ok)
fail_ss = sum(str(df.loc[t, f"passed_ss_{ic}"]).upper() != "TRUE" for t in q_ok)
expect(f"{['No', 'One', 'Two', 'Three', 'Four', 'Five', 'Six', 'Seven', 'Eight', 'Nine'][fail_ss]} of the {nominal} reach q",
       "q-passing interaction taxa that fail the sensitivity check")

# Genus-level forest plot: taxa named no further than genus.
gl = main.taxon.map(lambda t: [p.strip() for p in t.split(";") if p.strip() and not p.strip().endswith("__")])
_gp = gl.map(lambda ps: bool(ps) and not ps[-1].startswith('s__'))
expect(f"limited to the {int(_gp.sum())} taxa", "genus-level taxa")
_gq = [t for t in main.taxon[_gp.values] if qj.get(t, 1) < 0.05]
_words = ['No', 'One', 'Two', 'Three', 'Four', 'Five', 'Six', 'Seven', 'Eight', 'Nine']
expect(f"{_words[len(_gq)]} of the {_words[fail_ss].lower()} interactions with q &lt; 0.05 are among them", "genus-level interaction q count")
if sorted(t.rstrip("; _s").split("__")[-1] for t in _gq) != ["CAG-510", "Duncaniella", "Dysosmobacter"]:
    misses.append(f"genus-level q < 0.05 interactions differ from the text: {_gq}")
expect(f"; {_fail_ix} on the interaction itself", "finding interaction count")
expect(f"{nominal} of {len(df)} taxa have a nominal interaction where about {round(0.05 * len(df))} are expected by chance, "
       f"{_fail_ix} of them reach q &lt; 0.05 without passing", "synthesis interaction counts")
expect(f"{nominal} of {len(df)} taxa nominal where about {round(0.05 * len(df))} are expected by chance; {_fail_ix} reach q &lt; 0.05, none robust",
       "five-ways ANCOM-BC2 verdict")
expect(f"expected by chance ({_fail_ix} reach q &lt; 0.05, none robust)", "key finding interaction counts")

# ---- structural zeros: presence table and finding ----
import check_structural_zeros  # noqa: E402

sz_miss, szdf = check_structural_zeros.check(html, text)
misses += sz_miss
best = szdf[szdf.p == szdf.p.min()]
if set(best.taxon.map(lambda t: t.split(";")[-2])) != {"g__COE1"} or set(best.present_in) != {"IR"}:
    misses.append(f"strongest presence test is not COE1 in IR: {best[['taxon', 'present_in']].values.tolist()}")
expect(f"COE1 in {best.n_present.iloc[0]} of {best.n_group.iloc[0]} irradiated control-diet mice against none of the two "
       f"groups it is compared with, has p = {best.p.iloc[0]:.3f}", "COE1 presence")
for label, other in (("SCFAs+IR vs IR", "IR"), ("SCFAs+IR vs SCFAs", "unirradiated SCFAs mice")):
    d = szdf[szdf.contrast == label]
    expect(f"{int((d.absent_from == 'SCFAs+IR').sum())} of the {len(d)} against {other}", f"direction {label}")
expect(f"a taxon in all eight mice of one group and none of the other gives p = {fisher_exact([[8, 0], [0, 8]])[1]:.5f}",
       "design floor")

fig_miss, n_views = check_structural_zeros.check_figures(html, os.path.join(os.path.dirname(THIS_DIR), "comparison-Diet-x-Radiation-factorial"))
misses += fig_miss
if n_views == 0:
    misses.append("no view draws structural zeros")

# ---- presence interaction (presence_interaction.R) ----
def fmtp(p):
    return f"{p:.2g}" if p < 0.01 else (f"{p:.3f}" if p < 0.1 else f"{p:.2g}")


_pi = pd.read_csv(os.path.join(TAB, "ANCOMBC2_factorial_Greengenes2_presence_interaction.tsv"), sep="\t").set_index("term")
_ia, _pd = _pi.loc["Diet:Radiation"], _pi.loc["cell"]
expect(f"the presence interaction p = {fmtp(_ia.p)}", "presence interaction p")
expect(f"(R² {_ia.effect * 100:.0f}%)", "presence interaction R2")
expect(f"PERMDISP p = {fmtp(_pd.p)}", "presence PERMDISP p")
_nt = len(pd.read_csv(os.path.join(TAB, "ANCOMBC2_factorial_Greengenes2.tsv"), sep="\t"))
_nu = len(pd.read_csv(os.path.join(TAB, "ANCOMBC2_factorial_Greengenes2_untestable_absent_from_a_cell.tsv"), sep="\t"))
expect(f"The plot covers the {_nt} taxa seen in every group. The {_nu} absent", "forest scope")
expect(f"{_nt + _nu} taxa are seen in at least 4 of the 32 mice", "presence taxon set")
_bc = pd.read_csv(os.path.join(TAB, "ANCOMBC2_factorial_Greengenes2_presence_by_cell.tsv"), sep="\t", keep_default_na=False)
_bc = {(r[0], r[1]): r.taxa_present for r in _bc.itertuples(index=False)}
for (_a, _b), _c in [((('Control', 'None'), ('Control', 'IR')), None), ((('SCFAs', 'None'), ('SCFAs', 'IR')), None)]:
    expect(f"{_bc[_a]:.1f} to {_bc[_b]:.1f}" if _c is None else f"{_bc[_a]:.1f}" , f"presence count {_a}->{_b}")

_fr = pd.read_csv(os.path.join(TAB, "ANCOMBC2_factorial_Greengenes2_presence_firth.tsv"), sep="\t")
expect(f"{int((_fr.p_interaction < 0.05).sum())} of the {len(_fr)} taxa have raw p &lt; 0.05", "Firth counts")
expect(f"two opposite presence switches narrowly miss (q = {_fr.q_interaction.min():.2f})", "Firth table row")
expect(f"miss correction (q = {_fr.q_interaction.min():.2f})", "Firth q in findings")
expect(f"switch presence in opposite directions (q = {_fr.q_interaction.min():.2f})", "Firth q in synthesis")
expect(f"(p = {_fr.p_interaction.iloc[0]:.4f} and {_fr.p_interaction.iloc[1]:.4f})", "Firth two strongest p")
if not ("flavefaciens_B" in _fr.taxon.iloc[0] and _fr.taxon.iloc[1].endswith("g__COE1;s__")) or (_fr.q_interaction < 0.05).any():
    misses.append("Firth: the two strongest are not R. flavefaciens_B and genus COE1, or one passes correction")
_un = untestable.set_index("taxon"); _cc = [c for c in _un.columns if "(n=" in c]
_cell = lambda k: {c.split(" (")[0]: int(_un.loc[[t for t in _un.index if t.endswith(k)][0], c]) for c in _cc}
if _cell("flavefaciens_B") != {"Control": 4, "SCFAs": 0, "IR": 0, "SCFAs+IR": 4}:
    misses.append(f"R. flavefaciens_B presence differs from the text: {_cell('flavefaciens_B')}")
if _cell("g__COE1;s__") != {"Control": 0, "SCFAs": 3, "IR": 5, "SCFAs+IR": 0}:
    misses.append(f"genus COE1 presence differs from the text: {_cell('g__COE1;s__')}")
# Chance calibration of presence calls and complete separations (structural_zero_evidence.py).
_cal = pd.read_csv(os.path.join(TAB, f"{PRE}_structural_zeros_calibration.tsv"), sep="\t").set_index(["kind", "contrast"])
_a, _s = _cal.loc[("all presence calls", "all contrasts")], _cal.loc[("all presence calls", "SCFAs+IR vs SCFAs")]
_f = _cal.loc[("complete separations", "all contrasts")]
expect(f"there are {int(_a.observed)} presence calls where shuffled labels give {_a.expected:.1f} on average", "calibration all")
expect(f"({int(_s.observed)} against {_s.expected:.1f})", "calibration SCFAs+IR vs SCFAs")
expect(f"gives p = {_a.p_perm:.3f} (irradiation on SCFAs alone, {_s.p_perm:.3f})", "calibration p")
expect(f"({int(_f.observed)} complete separations, where chance gives {_f.expected:.3f})", "calibration complete")
expect(f"{int(_a.observed)} presence calls against {_a.expected:.1f} by chance", "summary calibration")
_irc = szdf[szdf.contrast == "IR vs Control"]
expect(f"{int((_irc.absent_from == 'Control').sum())} of the {len(_irc)} between irradiated and unirradiated control mice", "IR vs Control direction")
# ANCOM-BC2 finding: counts and the two taxa that move in opposite directions.
expect(f"{_fail['SCFAs vs Control']} between the diets in unirradiated mice", "finding diet count")
expect(f"and {_fail['SCFAs+IR vs IR']} in irradiated mice; {_fail['IR vs Control']} with irradiation on the control diet", "finding counts")
expect(f"and {_fail['SCFAs+IR vs SCFAs']} on the SCFAs diet", "finding SCFAs count")
_l2 = lambda suf, coef, k: pd.read_csv(os.path.join(TAB, f"{PRE}{suf}.tsv"), sep="\t").set_index("taxon").pipe(
    lambda t: t.loc[[x for x in t.index if x.endswith(k)][0], f"lfc_{coef}"] / LN2)
expect(f"falls {-_l2('', 'RadiationIR', 'Acetatifactor sp011959105'):.2f} log2 with irradiation", "Acetatifactor IR on control")
expect(f"is {_l2('_radrefIR', 'DietSCFAs', 'Acetatifactor sp011959105'):.2f} log2 higher in irradiated SCFAs", "Acetatifactor SCFAs+IR vs IR")
expect(f"Dysosmobacter rises {_l2('', 'RadiationIR', ';g__Dysosmobacter'):.2f} log2", "Dysosmobacter IR on control")
expect(f"is {-_l2('_radrefIR', 'DietSCFAs', ';g__Dysosmobacter'):.2f} log2 lower in irradiated SCFAs", "Dysosmobacter SCFAs+IR vs IR")
for _h in ("How taxa are tested", "Structural zeros: taxa present in some groups", "ANCOM-BC2: abundance in the 104 taxa", "Both together"):
    if _h not in html:
        misses.append(f"DA part heading missing: {_h}")
_order = [html.find(_h) for _h in ("How taxa are tested", "Structural zeros: taxa present", "ANCOM-BC2: abundance in", "Both together", "<h2>Summary</h2>")]
if _order != sorted(_order):
    misses.append("DA parts are out of order")

if _fr.p_interaction.isna().any():
    misses.append("Firth: some taxa have no interaction p-value")

# ---- PICRUSt2 (run_picrust2.sh), per-feature LinDA via ggpicrust2 at MetaCyc, KO and EC
# (run_ggpicrust2_linda.R, all features and --min-presence 8/12/16), PERMANOVA per level
# (picrust2_permanova.R) ----
import gzip  # noqa: E402
from statsmodels.stats.multitest import multipletests  # noqa: E402
_LV = ["MetaCyc", "KO", "EC"]; _CON = ["SCFAs vs Control", "SCFAs+IR vs IR", "IR vs Control", "SCFAs+IR vs SCFAs"]
_pm = pd.read_csv(os.path.join(TAB, "PICRUSt2_permanova.tsv"), sep="\t").set_index(["level", "term"])
_wn = pd.read_csv(os.path.join(OUT, "picrust2", "pipeline", "KO_metagenome_out", "weighted_nsti.tsv.gz"), sep="\t", index_col=0).weighted_NSTI
_nsti = pd.read_csv(os.path.join(OUT, "picrust2", "pipeline", "combined_marker_predicted_and_nsti.tsv.gz"), sep="\t", index_col=0).metadata_NSTI
if (_nsti > 2).any():
    misses.append("PICRUSt2: some ASVs exceed NSTI 2, but the text says none is dropped")
if len(_nsti) != 2573 or len(_wn) != 32:
    misses.append(f"PICRUSt2: {len(_nsti)} ASVs / {len(_wn)} mice, the text says 2,573 / 32")
expect(f"weighted NSTI runs from {_wn.min():.2f} to {_wn.max():.2f}", "NSTI range (caption)")

_L, _S = {}, {}
for _lv in _LV:
    for _mp, _suf in ((0, ""), (8, "_minpres8"), (12, "_minpres12"), (16, "_minpres16")):
        _d = pd.read_csv(os.path.join(TAB, f"PICRUSt2_{_lv}_linda{_suf}.tsv"), sep="\t", quoting=3)
        _L[(_lv, _mp)] = _d
        for _c in _CON:
            _x = _d[_d.contrast == _c]
            if not np.allclose(multipletests(_x.p, method="fdr_bh")[1], _x.q):
                misses.append(f"LinDA {_lv} {_mp} {_c}: q is not BH within the contrast")
            _sx = _x[_x.p < 0.05]
            _S[(_lv, _mp, _c)] = dict(n=len(_x), raw=len(_sx), down=int((_sx.log2fc < 0).sum()), qb=int((_x.q < 0.05).sum()),
                                     top=list(_x.sort_values("p").description.head(3)), med=_x.log2fc.median())
        if _mp == 16 and (_d.n_predicted < 16).any():
            misses.append(f"LinDA {_lv} >=16 fit contains features predicted in fewer than 16 mice")
        if _d.description.isna().any():
            misses.append(f"LinDA {_lv} {_mp}: features without a description")
_ex = lambda k: round(0.05 * k)
expect(f"{_S[('KO', 0, _CON[0])]['n']:,} KEGG orthologs (KO), {_S[('EC', 0, _CON[0])]['n']:,} enzymes (EC) and {_S[('MetaCyc', 0, _CON[0])]['n']} MetaCyc pathways", "feature counts")
# table rows
_pi = html.index("q &lt; 0.05 by filter: all / 8 / 12 / 16")
_rows = [[re.sub(r"<[^>]+>", " ", x).strip() for x in re.findall(r"<td[^>]*>(.*?)</td>", r)]
         for r in re.findall(r"<tr[^>]*>(.*?)</tr>", html[_pi:html.index("</tbody>", _pi)], re.S)]
_want = []
for _lv in _LV:
    for _c in _CON:
        a, w = _S[(_lv, 0, _c)], _S[(_lv, 16, _c)]
        qs = [_S[(_lv, mp, _c)]["qb"] for mp in (0, 8, 12, 16)]
        _want.append([_c, f"{a['raw']:,}", f"{w['raw']:,} / {w['n']:,} (~{_ex(w['n']):,})", " / ".join(f"{x:,}" for x in qs),
                      "filter-dependent" if any(qs) else "below threshold"])
_got = [r[-5:] for r in _rows]
if _got != _want:
    for g, w_ in zip(_got, _want):
        if g != w_:
            misses.append(f"LinDA table row: report {g}, tables {w_}")
    if len(_got) != len(_want):
        misses.append(f"LinDA table has {len(_got)} rows, expected {len(_want)}")
# PERMANOVA
_pr = lambda lv, t: _pm.loc[(lv, t)].effect * 100
_pv = lambda lv, t: _pm.loc[(lv, t)].p
_rng = lambda vals, fmt: f"{min(vals):{fmt}}–{max(vals):{fmt}}"
expect(f"diet explains {_rng([_pr(l, 'Diet') for l in _LV], '.1f')}% of the variation (p = {_rng([_pv(l, 'Diet') for l in _LV], '.3f')})", "synthesis diet")
expect(f"irradiation {_rng([_pr(l, 'Radiation') for l in _LV], '.1f')}% (p ≥ {min(_pv(l, 'Radiation') for l in _LV):.2f})", "synthesis irradiation")
expect(f"the interaction {_rng([_pr(l, 'Diet:Radiation') for l in _LV], '.1f')}% (p = {_rng([_pv(l, 'Diet:Radiation') for l in _LV], '.2f')})", "synthesis interaction")
expect(f"(PERMDISP p ≥ {min(_pv(l, 'group') for l in _LV):.2f})", "PERMDISP")
# findings
_w = {l: _S[(l, 16, "SCFAs+IR vs IR")] for l in _LV}; _v = {l: _S[(l, 16, "SCFAs+IR vs SCFAs")] for l in _LV}
_d0 = {l: _S[(l, 16, "SCFAs vs Control")] for l in _LV}; _i0 = {l: _S[(l, 16, "IR vs Control")] for l in _LV}
expect(f"SCFAs+IR against IR gives {_w['MetaCyc']['raw']} of {_w['MetaCyc']['n']} pathways, {_w['KO']['raw']:,} of {_w['KO']['n']:,} orthologs "
       f"and {_w['EC']['raw']} of {_w['EC']['n']:,} enzymes at raw p &lt; 0.05, where chance gives about "
       f"{_ex(_w['MetaCyc']['n'])}, {_ex(_w['KO']['n'])} and {_ex(_w['EC']['n'])}", "LinDA SCFAs+IR vs IR lean")
expect(f"{sum(_w[l]['down'] for l in _LV):,} of these {sum(_w[l]['raw'] for l in _LV):,} are lower in SCFAs+IR", "LinDA direction count")
expect(f"({_v['MetaCyc']['raw']}, {_v['KO']['raw']} and {_v['EC']['raw']})", "LinDA SCFAs+IR vs SCFAs counts")
expect(f"three ({_d0['MetaCyc']['raw']}, {_d0['KO']['raw']} and {_d0['EC']['raw']}; {_i0['MetaCyc']['raw']}, {_i0['KO']['raw']} and {_i0['EC']['raw']}, among features predicted in at least 16 mice)", "LinDA other contrasts")
if not all(_d0[l]["raw"] < _ex(_d0[l]["n"]) and _i0[l]["raw"] < _ex(_i0[l]["n"]) for l in _LV):
    misses.append("LinDA: a non-SCFAs+IR contrast is not below chance at >= 16")
expect(f"Among the strongest pathways are {_w['MetaCyc']['top'][0]}, {_w['MetaCyc']['top'][1]} and {_w['MetaCyc']['top'][2]}", "LinDA strongest pathways")
_q = {l: [_S[(l, mp, "SCFAs+IR vs IR")]["qb"] for mp in (0, 8, 12, 16)] for l in _LV}
expect(f"goes from {_q['MetaCyc'][0]} / {_q['KO'][0]} / {_q['EC'][0]} (pathways / orthologs / enzymes) with every feature to "
       f"{_q['MetaCyc'][3]} / {_q['KO'][3]:,} / {_q['EC'][3]} among features predicted in at least 16 mice", "LinDA q sweep")
expect(f"orthologs at 8 mice give {_S[('KO', 8, 'SCFAs+IR vs SCFAs')]['qb']:,} for SCFAs+IR against SCFAs and none in any other fit", "LinDA KO at 8")
if any(_S[("KO", mp, "SCFAs+IR vs SCFAs")]["qb"] for mp in (0, 12, 16)) or _q["KO"][1]:
    misses.append("LinDA KO: q counts at other filters are not zero as the text says")
_meds = [_S[(l, mp, c)]["med"] for l in _LV for mp in (0, 8, 12, 16) for c in ("SCFAs+IR vs IR", "SCFAs+IR vs SCFAs")]
expect(f"moves between {min(_meds):.2f} and {max(_meds):+.2f} log2", "LinDA median shift")
_e0 = _L[("EC", 0)]; _eh = _e0[(_e0.contrast == "SCFAs+IR vs IR") & (_e0.q < 0.05)]
if sorted(_eh.n_predicted) != [5, 32]:
    misses.append(f"EC q hits: predicted in {sorted(_eh.n_predicted)} mice, the text says 5 and 32")
_eh = _eh[_eh.n_predicted == 32].iloc[0]
expect(f"call on a well-predicted feature with every feature fitted is {_eh.description} (EC {_eh.feature}, predicted in all 32 mice), lower in SCFAs+IR ({_eh.log2fc:+.2f} log2, q = {_eh.q:.3f})", "EC well-predicted hit")
_k0 = _L[("KO", 0)]; _kh = _k0[(_k0.contrast == "SCFAs+IR vs IR") & (_k0.q < 0.05)]
if len(_kh) != 3 or set(_kh.n_predicted) != {5}:
    misses.append("KO q hits with every feature: not three KOs predicted in five mice")
_ks = _k0[(_k0.contrast == "SCFAs vs Control") & (_k0.p < 0.05)]
expect(f"{len(_ks)} orthologs differ between the diets without irradiation at raw p &lt; 0.05, {round(100 * (_ks.n_predicted < 16).mean())}% of them predicted in fewer than 16 mice", "KO sparse pattern")
expect(f"among well-predicted orthologs that contrast gives {_d0['KO']['raw']}, below chance", "KO sparse vs well-predicted")
_kab = pd.read_csv(os.path.join(OUT, "picrust2", "pipeline", "KO_metagenome_out", "pred_metagenome_unstrat.tsv.gz"), sep="\t", index_col=0)
_who = (_kab.loc[_ks[_ks.n_predicted < 16].feature] > 0).mean()
_nine = _who[_who > 0.5].index
_grp = pd.read_csv(os.path.join(OUT, "metadata_factorial.tsv"), sep="\t", index_col=0, keep_default_na=False).Group
if len(_nine) != 9 or _grp[_nine].value_counts().to_dict() != {"SCFAs": 4, "Control": 2, "IR": 2, "SCFAs+IR": 1}:
    misses.append(f"nine-mouse pattern: {_grp[_nine].value_counts().to_dict()}")
expect(f"({_rng([_pr(l, 'Diet') for l in _LV], '.0f')}%, p = {_rng([_pv(l, 'Diet') for l in _LV], '.3f')})", "summary diet")
expect(f"SCFAs+IR against IR gives {_w['MetaCyc']['raw']} pathways, {_w['KO']['raw']:,} KOs and {_w['EC']['raw']} enzymes at raw p &lt; 0.05 against about "
       f"{_ex(_w['MetaCyc']['n'])}, {_ex(_w['KO']['n'])} and {_ex(_w['EC']['n'])} by chance", "summary LinDA")
expect(f"{_rng([_pr(l, 'Diet:Radiation') for l in _LV], '.1f')}% of profile variation (p = {_rng([_pv(l, 'Diet:Radiation') for l in _LV], '.2f')})", "six-ways row")
if not (html.index("<h2>Differential abundance</h2>") < html.index("<h2>Predicted function (PICRUSt2)</h2>") < html.index("<h2>Summary</h2>")):
    misses.append("PICRUSt2 section is not between Differential abundance and Summary")
# Per-database subsections: PCA, heatmap and volcano of that database, a top-features table, a summary.
def _fmtp(x):
    if x < 0.001:
        return f"{x:.1e}"
    if x < 0.1:
        return f"{x:.4f}" if (x < 0.05) != (round(x, 3) < 0.05) else f"{x:.3f}"
    return f"{x:.2f}"
_LVN = {"MetaCyc": "MetaCyc pathways", "KO": "KEGG orthologs", "EC": "Enzymes"}
_vres = {}
for _cls, _lv in (("mc", "MetaCyc"), ("ko", "KO"), ("ec", "EC")):
    _m = re.search(rf'<input type="radio" name="vres"[^>]*value="{_cls}"', html)
    _i = html.find(f'<div class="view {_cls}">', html.index('name="vres"'))
    if not _m or _i < 0:
        misses.append(f"results panel: no {_lv} view"); continue
    _blk = html[_i:html.index("</tbody>", _i) + 8]
    _got = [[re.sub(r"<[^>]+>", "", x).strip() for x in re.findall(r"<td[^>]*>(.*?)</td>", r)]
            for r in re.findall(r"<tr[^>]*>(.*?)</tr>", _blk[_blk.index("<tbody>"):], re.S)]
    _want = []
    _all = _L[(_lv, 0)].set_index(["contrast", "feature"])
    for _c in _CON:
        _x = _L[(_lv, 16)]; _x = _x[_x.contrast == _c].sort_values("p").head(5)
        for _r in _x.itertuples():
            _qa = _fmtp(_all.loc[(_c, _r.feature)].q) if (_c, _r.feature) in _all.index else "—"
            _d = ("no name in either KO map" if _r.description == "not_found" else _r.description).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
            _want.append([_r.feature.replace("ko:", ""), _d, f"{_r.log2fc:+.2f}", _fmtp(_r.p), _fmtp(_r.q), _qa])
    _got6 = [g[-6:] for g in _got]
    if _got6 != _want:
        misses.append(f"{_lv} results view: table differs from the tables: {[(g, w) for g, w in zip(_got6, _want) if g != w][:2] or (len(_got6), len(_want))}")
    _wv, _vv = _S[(_lv, 16, "SCFAs+IR vs IR")], _S[(_lv, 16, "SCFAs+IR vs SCFAs")]
    _txt = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", _blk))
    for _w in (f"{_S[(_lv, 0, _CON[0])]['n']:,} predicted. Diet explains {_pm.loc[(_lv, 'Diet')].effect * 100:.1f}%",
               f"interaction {_pm.loc[(_lv, 'Diet:Radiation')].effect * 100:.1f}% (p = {_fmtp(_pm.loc[(_lv, 'Diet:Radiation')].p)})",
               f"LinDA gives {_wv['raw']:,} at raw p &lt; 0.05 for SCFAs+IR against IR ({_wv['down']:,} lower in SCFAs+IR; about {_ex(_wv['n']):,} expected by chance)",
               f"and {_vv['raw']:,} for SCFAs+IR against SCFAs; the other two contrasts give {_S[(_lv, 16, 'SCFAs vs Control')]['raw']:,} and {_S[(_lv, 16, 'IR vs Control')]['raw']:,}"):
        if _w not in _txt:
            misses.append(f"{_lv} results view: '{_w}' not found")
# Figure panels: each view holds its own database's image.
for _stem, _name in (("pca", "vpca"), ("heatmap", "vhm"), ("volcano", "vvol")):
    _p0 = html.index(f'name="{_name}"')
    for _cls, _lv in (("mc", "MetaCyc"), ("ko", "KO"), ("ec", "EC")):
        if not re.search(rf'<div class="view {_cls}">\s*<div class="one">\s*<img src="figures/picrust2_{_lv}_{_stem}\.png"', html[_p0:]):
            misses.append(f"panel {_stem}: no {_lv} view")
# Toggle wiring: each radio has its id and a label pointing at it.
for _name in ("vres", "vpca", "vhm", "vvol"):
    for _cls in ("mc", "ko", "ec"):
        if not re.search(rf'<input type="radio" name="{_name}" id="{_name}{_cls}" value="{_cls}"', html) or f'<label for="{_name}{_cls}">' not in html:
            misses.append(f"panel {_name}: {_cls} toggle not wired")
# Findings: exactly three, all synthesis across databases, after the panels.
_sec = html[html.index("<h2>Predicted function (PICRUSt2)</h2>"):html.index("</section>", html.index("<h2>Predicted function (PICRUSt2)</h2>"))]
_last_db_panel = max(_sec.index(f'name="{n}"') for n in ("vres", "vpca", "vhm", "vvol"))
if (len(re.findall(r'class="finding', _sec)) != 3
        or not _last_db_panel < _sec.index("What the three databases show together") < _sec.index("Acetate metabolism</h3>") < _sec.index('name="vac"')):
    misses.append("PICRUSt2 section: not database panels, then three synthesis findings, then the acetate panel")
# ---- acetate panel (picrust2_acetate.py) ----
_ac = pd.read_csv(os.path.join(TAB, "PICRUSt2_acetate_results.tsv"), sep="\t", quoting=3)
_acf = pd.read_csv(os.path.join(TAB, "PICRUSt2_acetate_features.tsv"), sep="\t", quoting=3)
_colmap = {"A": ("SCFAs vs Control", ""), "B": ("SCFAs+IR vs IR", ""), "C": ("IR vs Control", ""),
           "D": ("SCFAs+IR vs SCFAs", ""), "B′": ("SCFAs+IR vs IR", "_exclC2"), "D′": ("SCFAs+IR vs SCFAs", "_exclN3.2")}
for _r in _ac.itertuples():
    _sparse = _r.n_predicted < 16
    if _sparse != _r.role.startswith("Glyoxylate"):
        misses.append(f"acetate {_r.feature}: sparse flag and role disagree")
    _c, _suf = _colmap[_r.column]
    if _sparse and _suf:
        if not np.isnan(_r.log2fc):
            misses.append(f"acetate {_r.feature} {_r.column}: sparse feature has a one-mouse refit value")
        continue
    _src = _L[(_r.database, 0)] if _sparse else (_L[(_r.database, 16)] if not _suf else
            pd.read_csv(os.path.join(TAB, f"PICRUSt2_{_r.database}_linda_minpres16{_suf}.tsv"), sep="\t", quoting=3))
    _x = _src[(_src.feature == _r.feature) & (_src.contrast == _c)].iloc[0]
    if not (np.isclose(_x.log2fc, _r.log2fc) and np.isclose(_x.p, _r.p)):
        misses.append(f"acetate {_r.feature} {_r.column}: results file differs from LinDA table")
if sorted(zip(_acf.database, _acf.feature)) != sorted(set(zip(_ac.database, _ac.feature))):
    misses.append("acetate: results file and feature list differ")
_ACOL = ["A", "B", "C", "D", "B′", "D′"]
def _accell(r):
    if np.isnan(r.log2fc):
        return "—"
    return f"{r.log2fc:+.2f}p = {_fmtp(r.p)}"
for _cls, _db in (("mc", "MetaCyc"), ("ko", "KO"), ("ec", "EC")):
    if not re.search(rf'<input type="radio" name="vac" id="vac{_cls}" value="{_cls}"', html) or f'<label for="vac{_cls}">' not in html:
        misses.append(f"acetate panel: {_cls} toggle not wired")
    _i = html.find(f'<div class="view {_cls}">', html.index('name="vac"'))
    _blk = html[_i:html.index("</tbody>", _i)]
    if f'src="figures/picrust2_{_db}_acetate.png"' not in _blk:
        misses.append(f"acetate panel {_db}: figure missing from its view")
    _got = [[re.sub(r"<[^>]+>", "", x).strip() for x in re.findall(r"<td[^>]*>(.*?)</td>", r)]
            for r in re.findall(r"<tr>(.*?)</tr>", _blk[_blk.index("<tbody>"):], re.S)]
    _want = []
    _sub = _ac[_ac.database == _db]
    for (_role, _fid), _g in _sub.groupby(["role", "feature"], sort=False):
        _g = _g.set_index("column")
        _want.append([_g.label.iloc[0].replace("&", "&amp;"), str(int(_g.n_predicted.iloc[0]))] + [_accell(_g.loc[c]) for c in _ACOL])
    if _got != _want:
        misses.append(f"acetate panel {_db}: table differs: {[(g, w) for g, w in zip(_got, _want) if g != w][:1] or (len(_got), len(_want))}")
    for _c in re.findall(r'<td class="n">(.*?)</td>', _blk):
        _m = re.search(r"p = ([0-9.e-]+)", _c)
        if _m and (_c.startswith("<b>") != (float(_m.group(1)) < 0.05)):
            misses.append(f"acetate panel {_db}: bold does not match p < 0.05: {_c}")
_core = _ac[_ac.role.str.startswith("Acetate from acetyl-CoA")].dropna(subset=["log2fc"])
expect(f"change by at most {_core.log2fc.abs().max():.2f} log2 in any contrast or database (p ≥ {_core.p.min():.2f})", "acetate core route")
_g = lambda db, fid, col: _ac[(_ac.database == db) & (_ac.feature == fid) & (_ac.column == col)].iloc[0]
for _db, _fid, _txt in (("MetaCyc", "PWY-5100", "pyruvate fermentation to acetate and lactate ("), ("KO", "ko:K00198", "Wood–Ljungdahl CO dehydrogenase ("),
                        ("KO", "ko:K14138", "acetyl-CoA synthase ("), ("KO", "ko:K01938", "formate–tetrahydrofolate ligase ("),
                        ("KO", "ko:K01034", "acetate CoA-transferase ("), ("KO", "ko:K00158", "pyruvate oxidase ("), ("KO", "ko:K07034", "acetate uptake transporter (")):
    _x = _g(_db, _fid, "B")
    expect(f"{_txt}{_x.log2fc:+.2f}" + (" log2" if _fid == "PWY-5100" else "") + f", p = {_fmtp(_x.p)})", f"acetate {_fid}")
expect(f"hold without mouse C2 (p = {_fmtp(_g('MetaCyc', 'PWY-5100', 'B′').p)} and {_fmtp(_g('KO', 'ko:K01938', 'B′').p)})", "acetate without C2")
# PERMANOVA table, one row per database.
_pi = html.index("<th class=\"n\">Diet: R² (p)</th>")
_prow = [[re.sub(r"<[^>]+>", " ", x).strip() for x in re.findall(r"<td[^>]*>(.*?)</td>", r)]
         for r in re.findall(r"<tr[^>]*>(.*?)</tr>", html[_pi:html.index("</tbody>", _pi)], re.S)]
for _lv, _r in zip(_LV, _prow):
    _w = [f"{_pm.loc[(_lv, t)].effect * 100:.1f}% ({_fmtp(_pm.loc[(_lv, t)].p)})" for t in ("Diet", "Radiation", "Diet:Radiation")] + [_fmtp(_pm.loc[(_lv, "group")].p)]
    if _r[1:] != _w:
        misses.append(f"PERMANOVA table {_lv}: report {_r[1:]}, tables {_w}")
# Heatmap selection: 10 strongest per contrast at >= 16 mice, pooled.
_hm = pd.read_csv(os.path.join(TAB, "PICRUSt2_heatmap_features.tsv"), sep="\t", quoting=3)
for _lv in _LV:
    _d = _L[(_lv, 16)]
    _want = set().union(*[set(_d[_d.contrast == c].sort_values("p").feature.head(10)) for c in _CON])
    if set(_hm[_hm.level == _lv].feature) != _want:
        misses.append(f"heatmap {_lv}: feature set is not the 10 strongest per contrast")
    _i = html.index(f"figures/picrust2_{_lv}_heatmap.png")
    _cap = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", html[html.index("<figcaption>", _i):html.index("</figcaption>", _i)]))
    if f"pooled ({len(_want)} features)" not in _cap:
        misses.append(f"heatmap {_lv}: caption does not give {len(_want)} features")
# Heatmap captions name the dark SCFAs and IR columns: highest mean z over each level's heatmap rows.
_F = {"MetaCyc": "pathways_out/path_abun_unstrat.tsv.gz", "KO": "KO_metagenome_out/pred_metagenome_unstrat.tsv.gz",
      "EC": "EC_metagenome_out/pred_metagenome_unstrat.tsv.gz"}
for _lv in _LV:
    _ab = pd.read_csv(os.path.join(OUT, "picrust2", "pipeline", _F[_lv]), sep="\t", index_col=0)[_grp.index]
    _x = (_ab / _ab.sum()).loc[_hm[_hm.level == _lv].feature.unique()]
    _mz = _x.sub(_x.mean(1), axis=0).div(_x.std(1), axis=0).mean()
    if (_mz[_grp == "SCFAs"].idxmax(), _mz[_grp == "IR"].idxmax()) != ("N3.2", "C2"):
        misses.append(f"heatmap {_lv}: darkest SCFAs/IR columns are not N3.2/C2")
# One-mouse check: the highest comparison mouse and the refits without it.
_rel = _kab / _kab.sum()
def _top(c, g):
    f_ = _L[("KO", 16)]; f_ = f_[(f_.contrast == c) & (f_.p < 0.05)].feature
    vc = _rel.loc[f_, _grp.index[_grp == g]].idxmax(axis=1).value_counts()
    return vc.index[0], int(vc.iloc[0]), len(f_)
_ts, _ti = _top("SCFAs+IR vs SCFAs", "SCFAs"), _top("SCFAs+IR vs IR", "IR")
if (_ts[0], _ti[0]) != ("N3.2", "C2"):
    misses.append(f"one-mouse check: highest mice are {_ts[0]} / {_ti[0]}")
expect(f"IR mouse {_ti[0]} is the highest IR mouse on {_ti[1]:,} of the {_ti[2]:,} orthologs", "one-mouse C2 count")
expect(f"but that is SCFAs mouse {_ts[0]}, highest on {_ts[1]} of the {_ts[2]} orthologs", "one-mouse N3.2 count")
_X = {}
for _mo in ("N3.2", "C2"):
    for _lv in _LV:
        _d = pd.read_csv(os.path.join(TAB, f"PICRUSt2_{_lv}_linda_minpres16_excl{_mo}.tsv"), sep="\t", quoting=3)
        for _c in ("SCFAs+IR vs IR", "SCFAs+IR vs SCFAs"):
            _x = _d[_d.contrast == _c]
            _X[(_mo, _lv, _c)] = (len(_x), int((_x.p < 0.05).sum()), int((_x.q < 0.05).sum()))
_n3 = [_X[("N3.2", l, "SCFAs+IR vs SCFAs")] for l in _LV]
expect(f"without it the three give {_n3[0][1]}, {_n3[1][1]} and {_n3[2][1]}, below chance", "without N3.2")
if not all(x[1] < round(0.05 * x[0]) for x in _n3):
    misses.append("without N3.2: SCFAs+IR vs SCFAs is not below chance at every level")
_c2 = [_X[("C2", l, "SCFAs+IR vs IR")] for l in _LV]
expect(f"without it the three databases still give {_c2[0][1]}, {_c2[1][1]} and {_c2[2][1]}, about four times chance", "without C2 raw")
expect(f"and to {_c2[0][2]} / {_c2[1][2]} / {_c2[2][2]} without the dominant IR mouse", "without C2 q")
if not all(x[1] > 3 * round(0.05 * x[0]) for x in _c2):
    misses.append("without C2: SCFAs+IR vs IR is not several times chance at every level")

if misses:
    print("\n".join(misses))
    sys.exit(1)
print(f"report checks pass: {n_fig} figures, all quoted numbers found")
