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
for label, other in (("SCFAs+IR vs IR", "IR"), ("SCFAs+IR vs SCFAs", "unirradiated SCFAs mice")):
    d = szdf[szdf.contrast == label]
    expect(f"{int((d.absent_from == 'SCFAs+IR').sum())} of the {len(d)} against {other}", f"direction {label}")

fig_miss, n_views = check_structural_zeros.check_figures(html, os.path.join(os.path.dirname(THIS_DIR), "comparison-Diet-x-Radiation-factorial"))
misses += fig_miss
# The barplot/cladogram panels are shelved (shelved/both_together_barplot_cladogram.html);
# heat trees draw the structural zeros instead. Either kind of view counts.
ht_miss, n_ht = check_structural_zeros.check_heattrees(html, os.path.join(os.path.dirname(THIS_DIR), "comparison-Diet-x-Radiation-factorial"))
misses += ht_miss
n_views += n_ht
if n_ht and n_ht != 10:
    misses.append(f"{n_ht} heat-tree views, expected 10 (4 contrasts and the interaction, x species and genus)")
if n_views == 0:
    misses.append("no view draws structural zeros")

# ---- structural zeros: descriptive only (no presence statistics since 2026-10-01) ----
_nt = len(pd.read_csv(os.path.join(TAB, "ANCOMBC2_factorial_Greengenes2.tsv"), sep="\t"))
_nu = len(pd.read_csv(os.path.join(TAB, "ANCOMBC2_factorial_Greengenes2_untestable_absent_from_a_cell.tsv"), sep="\t"))
expect(f"The plot covers the {_nt} taxa seen in every group. The {_nu} absent", "forest scope")
expect(f"{_nt + _nu} taxa are seen in at least 4 of the 32 mice", "DA taxon set")
# Taxa per mouse by group: a count, written by presence_interaction.R beside its (retired) tests.
_bc = pd.read_csv(os.path.join(TAB, "ANCOMBC2_factorial_Greengenes2_presence_by_cell.tsv"), sep="\t", keep_default_na=False)
_bc = {(r[0], r[1]): r.taxa_present for r in _bc.itertuples(index=False)}
for _a, _b in ((("Control", "None"), ("Control", "IR")), (("SCFAs", "None"), ("SCFAs", "IR"))):
    expect(f"{_bc[_a]:.1f} to {_bc[_b]:.1f}", f"taxa per mouse {_a}->{_b}")
_un = untestable.set_index("taxon"); _cc = [c for c in _un.columns if "(n=" in c]
_cell = lambda k: {c.split(" (")[0]: int(_un.loc[[t for t in _un.index if t.endswith(k)][0], c]) for c in _cc}
if _cell("flavefaciens_B") != {"Control": 4, "SCFAs": 0, "IR": 0, "SCFAs+IR": 4}:
    misses.append(f"R. flavefaciens_B presence differs from the text: {_cell('flavefaciens_B')}")
if _cell("g__COE1;s__") != {"Control": 0, "SCFAs": 3, "IR": 5, "SCFAs+IR": 0}:
    misses.append(f"genus COE1 presence differs from the text: {_cell('g__COE1;s__')}")
if any((_un[_cc].values == 0).any(axis=1) & (_un[_cc].values == 8).any(axis=1)):
    misses.append("the text says no taxon is in every mouse of one group and none of another; one is")
expect("No taxon is in every mouse of one group and in none of the other", "no complete separation")
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


# ---- PICRUSt2 (run_picrust2.sh, picrust2_permanova.R, run_picrust2_maaslin3.R, plot_picrust2_targets.py):
# check_picrust2.py; the LinDA-era block is kept in the scratchpad as reza_check_report_linda.py ----
import check_picrust2  # noqa: E402
_pm_miss, _ko_nom, _ko_ch = check_picrust2.check(html)
misses += _pm_miss
if "LinDA" in html:
    misses.append("the report mentions LinDA, which it does not use")

# ---- ALDEx3 check on the interaction ----
_IC = "DietSCFAs:RadiationIR"
_anc = pd.read_csv(os.path.join(TAB, f"{PRE}.tsv"), sep="\t").set_index("taxon")
_al = {}
for _run in ("gamma0.5", "gamma0", "gamma0.5_seed7", "gamma0.5_ns2000"):
    _al[_run] = pd.read_csv(os.path.join(TAB, f"ALDEx3_factorial_Greengenes2_{_run}.tsv"), sep="\t").set_index("taxon")
    if set(_al[_run].index) != set(_anc.index):
        misses.append(f"ALDEx3 {_run}: taxa differ from the ANCOM-BC2 default fit")
_hits = lambda a: set(a.index[a[f"q_{_IC}"] < 0.05])
for _run in ("gamma0.5_seed7", "gamma0.5_ns2000"):
    if _hits(_al[_run]) != _hits(_al["gamma0.5"]):
        misses.append(f"ALDEx3 hit set changes in the {_run} rerun")

_x = pd.read_csv(os.path.join(TAB, "DA_interaction_intersection_Greengenes2.tsv"), sep="\t").set_index("taxon")
_a5 = _al["gamma0.5"].reindex(_x.index)
_lfc = _anc[f"lfc_{_IC}"].reindex(_x.index) / LN2
_both = (_anc[f"q_{_IC}"].reindex(_x.index) < 0.05) & (_a5[f"q_{_IC}"] < 0.05) & ((_lfc > 0) == (_a5[f"est_{_IC}"] > 0))
if not (np.allclose(_x.ancombc2_log2fc, _lfc) and np.allclose(_x["aldex3_log2fc_gamma0.5"], _a5[f"est_{_IC}"])
        and np.allclose(_x["aldex3_q_gamma0"], _al["gamma0"][f"q_{_IC}"].reindex(_x.index))
        and (_x.passes_both == _both).all() and (_x.n_zero_mice == _a5.n_zero_mice).all()):
    misses.append("intersection table is stale against the ANCOM-BC2 and ALDEx3 fits")

_fx = lambda v: f"{v:+.2f}".replace("-", "−")
_H = "Checked with ALDEx3: does the interaction hold"
_xh = _x[_x.ancombc2_q < 0.05].sort_values("ancombc2_q")
_tr = table_rows(_H)
if len(_tr) != len(_xh):
    misses.append(f"ALDEx3 table has {len(_tr)} rows, expected {len(_xh)}")
for _row, (_t, _r) in zip(_tr, _xh.iterrows()):
    _fq = lambda v: f"{v:.2f}" if v >= 0.1 else sig2(v)
    _want = [str(int(_r.n_zero_mice)), _fx(_r.ancombc2_log2fc), _fq(_r.ancombc2_q), _fx(_r["aldex3_log2fc_gamma0.5"]),
             _fq(_r["aldex3_q_gamma0.5"]), _fq(_r["aldex3_q_gamma0"]), "yes" if _r.passes_both else "no"]
    _name = _r.label.partition(" (")[0][3:]
    if not _row[0].startswith(_name) or _row[1:] != _want:
        misses.append(f"ALDEx3 table row for {_name}: {_row} != {_want}")
_body = html[html.index("<tbody>", html.index(_H)):html.index("</tbody>", html.index(_H))]
for _tr_html, (_t, _r) in zip(re.findall(r"<tr>(.*?)</tr>", _body, re.S), _xh.iterrows()):
    _tds = re.findall(r"<td[^>]*>(.*?)</td>", _tr_html)
    for _ci, _q in ((3, _r.ancombc2_q), (5, _r["aldex3_q_gamma0.5"]), (6, _r["aldex3_q_gamma0"])):
        if ("<b>" in _tds[_ci]) != (_q < 0.05):
            misses.append(f"ALDEx3 table bold does not match q < 0.05 for {_t} column {_ci}")

_zf = _x[_x.n_zero_mice == 0]
_off = abs((_zf["aldex3_log2fc_gamma0.5"] - _zf.ancombc2_log2fc).median())
_sp = _x[_x.n_zero_mice > 8]
_dun = _x[_x.label == "g__Duncaniella"].iloc[0]
_sph = _xh[_xh.n_zero_mice > 0]
_W = ['no', 'one', 'two', 'three', 'four', 'five', 'six', 'seven', 'eight']
expect(f"For the {len(_zf)} taxa present in all 32 mice, ALDEx3 ranks the interaction in exactly the order", "ALDEx3 zero-free count")
expect(f"every estimate {_off:.2f} log2 higher", "ALDEx3 offset")
expect(f"is {_fx(_dun.ancombc2_log2fc)} log2 in ANCOM-BC2 and {_fx(_dun['aldex3_log2fc_gamma0.5'])} in ALDEx3, "
       f"raw p = {sig2(_al['gamma0.5'].loc[_dun.name, f'p_{_IC}'])} but q = {sig2(_dun['aldex3_q_gamma0.5'])}", "ALDEx3 Duncaniella")
expect(f"The other {_W[len(_sph)]} hits have a zero count in {int(_sph.n_zero_mice.min())} to {int(_sph.n_zero_mice.max())} of the 32 mice",
       "ALDEx3 zero range of sparse hits")
expect(f"For {_W[int(_sph.sign_agrees.sum())]} of them ALDEx3 gives a similar estimate", "ALDEx3 same-sign count")
expect(f"across the {len(_sp)} taxa with a zero in more than 8 mice, the median standard error is "
       f"{_al['gamma0.5'].loc[_sp.index, f'se_{_IC}'].median():.2f} log2 in ALDEx3 against "
       f"{(_anc.loc[_sp.index, f'se_{_IC}'] / LN2).median():.2f} in ANCOM-BC2", "ALDEx3 SE medians")
_flip = [("genus-level " if _l.startswith("g__") else "") + _l.partition(" (")[0][3:] for _l in _sph.label[~_sph.sign_agrees]]
_tight = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", "", html))
if f"For the other two, {' and '.join(_flip)}, the sign changes" not in _tight:
    misses.append(f"ALDEx3 sign flips: {_flip} not named as in the table")
_q5, _q0 = sig2(_x["aldex3_q_gamma0.5"].min()), sig2(_x["aldex3_q_gamma0"].min())
expect(f"(smallest q {_q5} and {_q0})", "ALDEx3 smallest q, both gammas")
expect(f"ALDEx3 puts no interaction below q = 0.05 (smallest q {_q5})", "ALDEx3 summary row")
expect(f"No taxon at q &lt; 0.05 (smallest q {_q5}); agrees with ANCOM-BC2 on the {len(_zf)} taxa in every mouse", "ALDEx3 seven-ways row")
expect(f"{len(_xh)} in ANCOM-BC2, none in ALDEx3", "ALDEx3 figure caption counts")
if int(_x.passes_both.sum()) != 0 or int((_x["aldex3_q_gamma0.5"] < 0.05).sum()) != 0:
    misses.append("the text says no taxon passes in ALDEx3, but the table has some")
# ---- presence dot map (plot_presence_dots_factorial.py): first panel of Both together ----
_un = pd.read_csv(os.path.join(TAB, f"{PRE}_untestable_absent_from_a_cell.tsv"), sep="\t")
_ntest = len(pd.read_csv(os.path.join(TAB, f"{PRE}.tsv"), sep="\t"))
_cells = {g: next(c for c in _un.columns if c.startswith(f"{g} (n=")) for g in ("Control", "IR", "SCFAs", "SCFAs+IR")}
for _gp in (False, True):
    _stem = f"{PRE}{'_gplus' if _gp else ''}_presence_dots"
    _sd = pd.read_csv(os.path.join(OUT, "figures", f"{_stem}.tsv"), sep="\t").iloc[0]
    _u = _un[[not _t.split(";")[-1].strip()[3:] or not _t.split(";")[-1].strip().startswith("s__") for _t in _un.taxon]] if _gp else _un
    _want = {f"absent_{g}": int((_u[c] == 0).sum()) for g, c in _cells.items()}
    _want["structural_zeros"] = len(_u)
    if not _gp:
        _want["taxa"] = _ntest + len(_un)
    for _k, _v in _want.items():
        if int(_sd[_k]) != _v:
            misses.append(f"{_stem}: sidecar {_k} = {_sd[_k]}, recounted {_v}")
    if f'src="figures/{_stem}.png"' not in html:
        misses.append(f"{_stem}: not in the report")
        continue
    _cap = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", html[html.index(f'src="figures/{_stem}.png"'):].split("</figcaption>")[0]))
    if _gp:
        _w = [f"for the {int(_sd.taxa)} taxa named no further than genus; {int(_sd.structural_zeros)} of them are structural zeros"]
    else:
        _gaps = sorted(((int(_sd[f"absent_{g}"]), g) for g in _cells), reverse=True)
        _w = [f"the {_ntest} tested and the {int(_sd.structural_zeros)} structural zeros: {int(_sd.taxa)} in all",
              "no dot in that column: " + ", ".join(f"{n} in {g}" for n, g in _gaps[:-1]) + f" and {_gaps[-1][0]} in {_gaps[-1][1]}"]
    for _x in _w:
        if _x not in _cap:
            misses.append(f"{_stem}: caption lacks '{_x}'")
if ">Both together</h3>" in html:
    misses.append("the DA section has no 'Both together' part any more")
if "_untestable_heatmap.png" in html:
    misses.append("the per-mouse structural-zero heatmap is still in the report")

# ---- community composition: group means of the species that differ most, shares only ----
_c = counts.loc[counts.index.intersection(tax.index)]
_paths = tax.loc[_c.index].map(lambda t: ";".join([p for p in (x.strip() for x in t.split(";")) if p][:7]))
_c = _c.groupby(_paths).sum()
_meta = pd.read_csv(os.path.join(os.path.dirname(OUT), "metadata.tsv"), sep="\t", index_col=0)
_c = _c[_meta.index.intersection(_c.columns)]
_rel = _c / _c.sum() * 100
_G = ["Control", "SCFAs", "IR", "SCFAs+IR"]
_mean = pd.DataFrame({g: _rel.loc[:, _meta.loc[_rel.columns, "Group"] == g].mean(axis=1) for g in _G})
_mean = _mean[~_mean.index.str.endswith("g__;s__")]
_top = (_mean.max(axis=1) - _mean.min(axis=1)).sort_values(ascending=False).index[:12]
_cs = html[html.index("<h2>Community composition</h2>"):html.index("</section>", html.index("<h2>Community composition</h2>"))]
_rows = re.findall(r"<tr><td>.*?</tr>", _cs)
if len(_rows) != 12:
    misses.append(f"composition table has {len(_rows)} rows, expected 12")
for _r, _t in zip(_rows, _top):
    _want = [("&lt;0.1%" if v < 0.05 else f"{v:.1f}%") for v in _mean.loc[_t, _G]]
    _got = re.findall(r'<td class="n">(.*?)</td>', _r)
    if _got != _want:
        misses.append(f"composition row {_t.split(';')[-1]}: report {_got}, tables {_want}")
_ext = sum(_mean.loc[t, _G].idxmax() == "SCFAs+IR" or _mean.loc[t, _G].idxmin() == "SCFAs+IR" for t in _top)
expect(f"Irradiated SCFAs mice are the extreme group for {_ext} of the 12", "composition extremes")
expect(f"extreme group for {_ext} of the 12 species that differ most", "summary composition extremes")
if re.search(r"\bq\b|\bq =|raw p", re.sub(r"<[^>]+>", " ", _cs)):
    misses.append("Community composition quotes a test result; it precedes the DA methods")

# Composition bars: one toggled view per rank, Kingdom to Species, Genus shown by default.
for _rk in ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"):
    _v = f"rk-{_rk[0].lower()}"
    if (f'<div class="view {_v}">' not in _cs or f'taxa_factorial_Greengenes2_{_rk}_groups.png' not in _cs
            or f'value="{_v}"' not in _cs or f'.view.{_v}{{display:block}}' not in html):
        misses.append(f"composition rank panel: {_rk} view missing or not wired")
if 'value="rk-g" checked' not in _cs:
    misses.append("composition rank panel: Genus is not the default view")

# ---- Summary key findings: each number is one the sections above quote ----
_kf = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", html[html.index("Key findings</h3>"):]))
_body = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", html[:html.index("<h2>Summary</h2>")]))
for _n in ("q = 0.17", "q = 0.076", "27 of 104", f"{_ko_nom:,} KOs at raw p &lt; 0.05, where chance would give about {_ko_ch}",
           "56.4 against 64.0", "69.5 to 74.1", "3.0%"):
    if _n not in _kf:
        misses.append(f"key findings: '{_n}' missing")
for _n in ("q = 0.17", "q = 0.076", "27 of 104", "64.0 to 56.4", "69.5 to 74.1", "3.0% (2.0–3.4%)"):
    if _n not in _body:
        misses.append(f"key findings quote '{_n}', which no section above states")
expect(f"irradiated SCFAs mice have the highest or lowest share for {_ext}.", "key finding composition extremes")
expect("gives the most taxa at q &lt; 0.05 (33)", "key finding SCFAs+IR calls")

# ---- section order: each analysis runs explanation, table, main finding, plot, caveats ----
def _order(section, steps):
    _s0 = html.index(section)
    _sec = html[_s0:html.index("</section>", _s0)]
    _p = [_sec.find(t) for t in steps]
    if -1 in _p or _p != sorted(_p):
        misses.append(f"{section}: out of order or missing: {[t for t, q in zip(steps, _p) if q == -1] or 'order'}")


_order("<h2>Alpha diversity</h2>", ["Two-way ANOVA on", "<thead>", "None survives multiple-testing correction",
                                    "diversity_alpha.png", "Three limits apply", 'class="finding'])
_order("<h2>Beta diversity</h2>", ["PCoA: how far apart", "A principal coordinates analysis", "the four groups overlap heavily",
                                   "diversity_beta_pcoa.png", "A PCoA is a picture", "PERMANOVA: do the groups differ",
                                   "PERMANOVA is <code>adonis2", "<th>Term</th>", "No term passes on any metric",
                                   "diversity_permanova_r2.png", "PERMANOVA has one blind spot", "PERMDISP: do the groups differ",
                                   "PERMDISP (vegan's", "<th>Grouping</th>", "Spread is equal in every row",
                                   "diversity_dispersion.png", 'class="finding'])
_order("<h2>Community composition</h2>", ["Mean relative abundance in each group", "<th>Species · Greengenes2</th>",
                                          "The same lineages dominate", "taxa_factorial_Greengenes2_Genus_groups.png",
                                          "A group mean can be carried"])
_order("<h2>Differential abundance</h2>", ["104 of them are seen", "Structural zeros: taxa present", "cannot enter a model of log",
                                           "Absent from first / second", "31 taxa are seen in at least 4 mice",
                                           "_presence_dots.png", "The share column is how much",
                                           "ANCOM-BC2: abundance in the 104 taxa", "Version 2.14.0", "Largest up / down",
                                           "Taxa reach q &lt; 0.05 in every row", "_heattree_DietSCFAs.png",
                                           "_interaction_forest.png", "_heattree_interaction.png", "ANCOM-BC2 corrects for differences",
                                           "Checked with ALDEx3", "Version 1.3.1", "Mice with a zero", "For the 17 taxa present in all 32 mice",
                                           "DA_interaction_concordance_Greengenes2.png", 'class="finding'])

# ---- every section: at most four findings, and they close it ----
# Sections are the <section> elements; each must open with its <h2>. The Summary alone may add up
# to two red findings ("finding no": easily misread results) after its four green or orange ones,
# with only a heading between the two groups (user, 2026-10-01).
for _m in re.finditer(r"<section[^>]*>(.*?)</section>", html, re.S):
    _sec = _m.group(1)
    _h = re.search(r"<h2>(.*?)</h2>", _sec)
    _name = _h.group(1) if _h else "(section without a heading)"
    if not _sec.strip().startswith("<h2>"):
        misses.append(f"{_name}: something precedes the section heading: "
                      f"{re.sub(r'<[^>]+>', ' ', _sec.strip()[:80])!r}")
    _boxes = list(re.finditer(r'<div class="finding([^"]*)">.*?</div>', _sec, re.S))
    _red = [b for b in _boxes if " no" in b.group(1)] if _name == "Summary" else []
    _other = [b for b in _boxes if b not in _red]
    if len(_other) > 4:
        misses.append(f"{_name}: {len(_other)} findings, at most 4")
    if len(_red) > 2:
        misses.append(f"{_name}: {len(_red)} red findings, at most 2")
    if _red and _other and _red[0].start() < _other[-1].start():
        misses.append(f"{_name}: a red finding comes before a green/orange one")
    if _boxes:
        _after = re.sub(r"\s+", " ", re.sub(r"<script.*?</script>|<[^>]+>", " ", _sec[_boxes[-1].end():], flags=re.S)).strip()
        if _after:
            misses.append(f"{_name}: content after its findings: {_after[:60]!r}")
        _rest = re.sub(r'<script.*?</script>|<div class="finding[^"]*">.*?</div>', "", _sec[_boxes[0].start():], flags=re.S)
        if _red and _other:   # a heading may introduce the red group
            _rest = re.sub(r"<h3[^>]*>.*?</h3>", "", _rest, flags=re.S)
        _between = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", _rest)).strip()
        if _between:
            misses.append(f"{_name}: findings are not together at the end ({_between[:60]!r} between them)")

if misses:
    print("\n".join(misses))
    sys.exit(1)
print(f"report checks pass: {n_fig} figures, all quoted numbers found")
